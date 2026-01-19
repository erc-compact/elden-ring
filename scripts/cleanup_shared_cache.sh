#!/bin/bash
# Cleanup unreferenced files in shared_cache
# This script scans the shared_cache directory and identifies files that are not
# referenced by any symlinks in runID-specific directories

set -euo pipefail

BASEDIR="${1:-/bscratch/fazal/ELDEN}"
SHARED_CACHE="${BASEDIR}/shared_cache"
DRY_RUN="${2:-true}"  # Default to dry run

if [[ ! -d "${SHARED_CACHE}" ]]; then
    echo "ERROR: Shared cache directory not found: ${SHARED_CACHE}"
    exit 1
fi

echo "================================================================"
echo "Shared Cache Cleanup Utility"
echo "================================================================"
echo "Base directory: ${BASEDIR}"
echo "Shared cache:   ${SHARED_CACHE}"
echo "Mode:           $([ "$DRY_RUN" = "true" ] && echo "DRY RUN (no files will be deleted)" || echo "LIVE (files will be deleted)")"
echo "================================================================"
echo ""

# Count statistics
total_files=0
orphaned_files=0
referenced_files=0
total_size=0
orphaned_size=0

echo "Scanning for orphaned cache files..."
echo ""

# Find all files in shared_cache
while IFS= read -r -d '' cache_file; do
    ((total_files++))

    # Get file size
    file_size=$(stat -f%z "$cache_file" 2>/dev/null || stat -c%s "$cache_file" 2>/dev/null || echo 0)
    ((total_size+=file_size))

    # Check if any symlink points to this file
    # Search in all runID directories (exclude shared_cache itself)
    refs=$(find "${BASEDIR}" -maxdepth 4 -type l \( -lname "*${cache_file}*" -o -lname "*$(basename ${cache_file})*" \) 2>/dev/null | grep -v shared_cache | wc -l)

    if [[ $refs -eq 0 ]]; then
        ((orphaned_files++))
        ((orphaned_size+=file_size))

        # Human-readable size
        if command -v numfmt &> /dev/null; then
            hr_size=$(numfmt --to=iec-i --suffix=B $file_size)
        else
            hr_size="${file_size} bytes"
        fi

        echo "ORPHANED [$hr_size]: ${cache_file}"

        # Delete if not in dry run mode
        if [[ "$DRY_RUN" != "true" ]]; then
            rm -f "${cache_file}"
            echo "  → DELETED"
        fi
    else
        ((referenced_files++))
    fi
done < <(find "${SHARED_CACHE}" -type f -print0)

echo ""
echo "================================================================"
echo "Scan Complete - Summary"
echo "================================================================"
echo "Total files scanned:     ${total_files}"
echo "Referenced files:        ${referenced_files}"
echo "Orphaned files:          ${orphaned_files}"

# Human-readable sizes
if command -v numfmt &> /dev/null; then
    echo "Total cache size:        $(numfmt --to=iec-i --suffix=B $total_size)"
    echo "Orphaned file size:      $(numfmt --to=iec-i --suffix=B $orphaned_size)"
else
    echo "Total cache size:        ${total_size} bytes"
    echo "Orphaned file size:      ${orphaned_size} bytes"
fi

echo "================================================================"

if [[ "$DRY_RUN" = "true" ]]; then
    echo ""
    echo "This was a DRY RUN - no files were deleted."
    echo "To actually delete orphaned files, run:"
    echo "  bash $0 ${BASEDIR} false"
else
    echo ""
    echo "Cleanup complete. ${orphaned_files} orphaned files were deleted."
fi

exit 0
