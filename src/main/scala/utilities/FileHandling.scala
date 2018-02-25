package utilities

import java.io.File
import java.util.Calendar
import scala.io.Source

/**
  * Author: Alex N. Salazar
  * Created on 14-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description: set of methods/traits for dealing with files (opening, iterators, parsing, etc.)
  */
object FileHandling {


  /**Method to get current time stamp*/
  def timeStamp = " (" + Calendar.getInstance().getTime() + ")"

  /**Method to get name of file*/
  def getFileName: File => String = file => file.getName.substring(0, file.getName.lastIndexOf("."))

  /**Method to open file with Iterator*/
  def openFileWithIterator(x: File, remove_header_comment_lines: Boolean = false): Iterator[String] = {
    if(remove_header_comment_lines) Source.fromFile(x).getLines.dropWhile(_.startsWith("#"))
    else Source.fromFile(x).getLines
  }

  /**Check whether directory exists and if it is valid. If it does not exist, create it.*/
  def checkOutputDirectory(inputDir: File) {
    if (inputDir.exists()) {
      assert(inputDir.isDirectory(), "The following directory does not exist or it is invalid: " + inputDir.getAbsolutePath())
    } else {
      inputDir.mkdir()
      assume(inputDir.exists(), "Could not create the following directory: " + inputDir.getAbsolutePath())
    }
  }

  /**Create directory and verify that it now exists*/
  def createDirectory(inputDir: File, message: String = "Could not create the following directory") {
    inputDir.mkdir()
    verifyDirectory(inputDir, message)
  }

  /**Check whether directory exists and if it is valid. If it does not exist, create it.*/
  def verifyDirectory(inputDir: File, message: String = "The following directory does not exist or it is invalid") {
    assert(inputDir.exists() && inputDir.isDirectory(), message + ": " + inputDir.getAbsolutePath())
  }

  /**Check whether directory exists and if it is valid. If it does not exist, create it.*/
  def verifyFile(inputFile: File, message: String = "The following directory does not exist or it is invalid") {
    assert(inputFile.exists() && inputFile.isFile(), message + ": " + inputFile.getAbsolutePath())
  }

  /**Method to clean a given directory*/
  def cleanDirectory(coverage_dir: File) {
    //make directory if it already exists
    println("Deleting files from the following subdirectory: " + coverage_dir.getAbsolutePath + " " + timeStamp)
    coverage_dir.listFiles().foreach(file => file.delete())
  }

  /**
    * Method to make SLURM headers
    * @param mem Allocated memory in Mb
    * @param time Allocated run time in hours
    * @param partition Requested partition
    * @param qos Requested type
    * @param output Full path of SLURM log output file
    * @param cpus Optional. Number of cpus available. Default is 1.
    */
  def makeSLURMheader(mem: Int, time: Int, partition: String, qos: String, output: String,
                      cpus: Int = 1): String = {
    //determine if additional SBATCH lines are needed if requesting multiple cpus
    val cpu_lines = {
      cpus match {
        //do nothing
        case 1 => ""
        //add cpus and number of tasks
        case _ => Seq(("#SBATCH --cpus-per-task=" + cpus), "#SBATCH --ntasks=1").mkString("\n")
      }
    }
    //create SLURM parameters and add 'set -x' flag
    Seq("#!/bin/sh",
      "#SBATCH --partition=" + partition + " --qos=" + qos,
      "#SBATCH --mem=" + mem,
      "#SBATCH --time " + time + ":00:00",
      "#SBATCH --out=" + output,
      cpu_lines).mkString("\n") + "\nset -x\n"
  }

  trait tLines{
    /**##################Following code below copied from Thomas Abeel ATK library###########################*/
    implicit def toLeft[String, File](left: String): Either[String, File] = Left(left)

    implicit def toRight[String, File](right: File): Either[String, File] = Right(right)

    def tLines(file: Either[String, File], skipComments: Boolean = true, skipBlank: Boolean = true): List[String] = {
      val ff = file.fold(new File(_), identity)
      if (!ff.exists()) {
        var parent = ff
        while (parent != null && !parent.exists()) {
          System.err.println("Invalid path: " + parent)
          parent = parent.getParentFile()
        }
      }
      Source.fromFile(ff).getLines().filterNot(x => x.startsWith("#") || x.isEmpty).toList
    }
    /**#####################################################################################################*/
  }

}
