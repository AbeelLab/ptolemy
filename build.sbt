name := "Ptolemy"

version := "0.1"

scalaVersion := "2.12.4"

libraryDependencies ++= Seq(
  "org.scalatest" %% "scalatest" % "3.0.4" % "test",
  "com.github.scopt" % "scopt_2.12" % "3.7.0",
  "com.sksamuel.avro4s" %% "avro4s-core" % "1.9.0",
  "com.github.jbellis" % "jamm" % "0.3.0",
  "io.verizon.quiver" % "core_2.12" % "7.0.18"
)