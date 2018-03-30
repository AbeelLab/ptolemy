package utilities

import java.io.{File, PrintWriter}

import scala.collection.immutable.HashMap
import scala.sys.process.ProcessLogger
import scala.sys.process._
import utilities.FileHandling.openFileWithIterator

/**
  * Author: Alex N. Salazar
  * Created on 14-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
trait MinimapUtils extends ORFalignments{

  /**
    * Method to peform pariwise ORF alignments and retain best alignments
    *
    * @param query
    * @param ref
    * @param gamma
    * @param hashmap_H
    * @return
    */
  def collectBestAlignments(ref: File, query: File,
                            gamma: Double,
                            kmerSize: Int,
                            minimizerWindow: Int,
                            hashmap_H: HashMap[Int,Set[Int]],
                            self: Boolean = false): HashMap[Int, Set[Int]] = {

    /**
      * Function to determine if the similarity of two ORF alignments are similar in sequence identity enough by
      * requiring the ratio of the #matches and length each ORF sequence to be >= sigma.
      *
      * @return Boolean
      */
    //def isSimilar: Array[String] => Boolean = alignment => computeSigma(alignment) >= sigma

    /**
      * Function to determine if the coverage of two ORF alignments are large enough by requiring the ratio of the
      * matches and length each ORF sequence to be >= gamma.
      *
      * @return Boolean
      */
    def isCovered: Array[String] => Boolean = alignment => computeGamma(alignment) >= gamma
    //command for running minimap2
    val command ={
      if(self) {
        Seq("minimap2", "-X", "-k", kmerSize.toString, "-w", minimizerWindow.toString, ref.getAbsolutePath,
          query.getAbsolutePath)
      }
      else
        Seq("minimap2", "-k", kmerSize.toString, "-w", minimizerWindow.toString, ref.getAbsolutePath, query.getAbsolutePath)
    }
    //capture stdout and stderr; note: using mutable stringbuilder, dirty but gets job done; may be optimized later
    var out = new StringBuilder
    var err = new StringBuilder
    val logger = ProcessLogger((o: String) => out.append(o+"\n"), (e: String) => err.append(e))
    //run command
    command ! logger
    //get alignments, and retain those that pass the sigma and gammas filters, add remaining to hashmap H
    val hashmap = out.toString.split("\n").foldLeft(hashmap_H) { case (map, line) => {
      //split line
      val alignment = line.split("\t")
      //alignment does not meet sigma and gamma quality thresholds
      if (!isCovered(alignment)) map
      else {
        //get current values for query ORF
        val current = map.getOrElse(alignment.head.toInt, Set[Int]())
        //update hashmap_H
        map + (alignment.head.toInt -> (current + alignment(5).toInt))
      }
    }
    }
    //sanity check
    assume(!hashmap.isEmpty, "Hashmap H for " + query.getName + ", ")
    //return hashmap
    hashmap
  }


}
