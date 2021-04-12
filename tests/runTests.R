if(ChemmineOB:::.supportedPlatform()){
  BiocGenerics:::testPackage("ChemmineOB")
}else{
  message("Not running tests on unsupported platform ", .Platform$OS, " ", .Platform$r_arch)
}
