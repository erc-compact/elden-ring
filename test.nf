process foo {
  input:
  val x
  output:
  tuple path('*fq'),path('*.ar'), path('*.cands')

  script:
  """
  #!/bin/bash
    echo "Hello World" > file.fq
    echo "cands" > file.cands
    touch file.ar
  """
}

process bar {
  debug true   
  input:
  path (file)
  script:
  """
  cat *.cands | head -n 50
  """
}

workflow {
    Channel.of('a', 'b', 'c', '1', '2', '3', '4', '5') \
    | foo \
    | collectFile(name: 'finished.txt', newLine : true) \
    | bar
}