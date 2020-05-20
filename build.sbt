name := "Ptolemy"

version := "0.1"

scalaVersion := "2.12.4"
scalacOptions += "-target:jvm-1.8"
javacOptions ++= Seq("-source", "1.8", "-target", "1.8")

libraryDependencies ++= Seq(
  "org.scalatest" %% "scalatest" % "3.0.4" % "test",
  "com.github.scopt" % "scopt_2.12" % "3.7.0",
  "com.sksamuel.avro4s" %% "avro4s-core" % "1.9.0",
  "io.verizon.quiver" % "core_2.12" % "7.0.18",
  "org.scala-lang.modules" %% "scala-parser-combinators" % "1.1.2",
)

// set the main class for packaging the main jar
mainClass in (Compile, packageBin) := Some("cli.Ptolemy")
mainClass in assembly := Some("cli.Ptolemy")

// set the main class for the main 'sbt run' task
mainClass in (Compile, run) := Some("cli.Ptolemy")
