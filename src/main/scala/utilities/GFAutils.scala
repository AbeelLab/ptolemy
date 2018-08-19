package utilities

import java.io.File

import quiver.{Context, Graph, LEdge, LNode, empty}
import utilities.FileHandling.openFileWithIterator
import scala.collection.immutable.HashSet

import scala.collection.immutable.{HashMap, HashSet}

/**
  * Author: Alex N. Salazar
  * Created on 16-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GFAutils {

  trait ConstructGFA {

    val getGFAHeader = "H\tHLGG"

    /**
      * Function to construct a Segment line for a given node
      *
      * @return String
      */
    def constructSegmentLine: Int => String = node => Seq("S", node, "N").mkString("\t")

    /**
      * Function to construct a Link line for a given pair of nodes
      *
      * @return
      */
    def constructLinkLine: (Int, Int) => String = (node, edge) =>
      Seq("L", node, "+", edge, "+", "1M").mkString("\t")

    /**
      * Function to create optional column of the number of genomes in node
      *
      * @return
      */
    def addGenomeCountField: Option[Int] => String = genomes =>
      "\tFC:i:" + {
        if (genomes == None) 1 else genomes.get
      }

  }

  trait ReadGFA {

    /**
      * Get the graph, the paths, and the genomes from a GFA file
      * @param gfa GFA file
      * @return
      */
    def loadGFA(gfa: File): (Map[Int, Set[Int]], Map[String, Seq[Int]], Map[String, Set[String]]) = {
      //iterate through each line and process only path or genome lines
      openFileWithIterator(gfa).foldLeft((Map[Int, Set[Int]](), Map[String, Seq[Int]](), Map[String, Set[String]]())) {
        case ((graph, paths, genomes), line) => {
          if (line.startsWith("S")) {
            val node = parseSegmentLine(line)
            (graph + (node -> graph.getOrElse(node, Set[Int]())), paths, genomes)
          }
          else if (line.startsWith("L")) {
            val (node, edge) = parseLinkLine(line)
            val current = graph.getOrElse(node, Set[Int]())
            (graph + (node -> (current + edge)), paths, genomes)
          }
          else if (line.startsWith("G")) {
            //split line
            val split = line.split("\t")
            assert(genomes.get(split(1)) == None)
            (graph, paths, genomes + (split(1) -> split(2).split(",").toSet))
          }
          else if (line.startsWith("P")) {
            //split line
            val split = line.split("\t")
            //assure genome has never been seen before
            assert(paths.get(split.head) == None)
            //add path to map
            (graph, paths + (split(1) -> split(2).split("""\W+""").toSeq.map(_.toInt)), genomes)
          } else (graph, paths, genomes)
        }
      }
    }

    /**
      * Function to load only the nodes in the GFA
      * @param gfa GFA file
      * @return List[Int]
      */
    def loadNodesGFA(gfa: File): List[Int] = {
      //iterate through each line and process only segment lines
      openFileWithIterator(gfa).foldLeft(List[Int]())((nodes, line) => {
        if(!line.startsWith("S")) nodes else nodes.:+(parseSegmentLine(line))
      })
    }

    /**
      * Get only the graph and the paths in a GFA file.
      * @param gfa GFA file
      * @return Map as ID -> Seq[Node IDs]

    def getSequencesPaths(gfa: File): (Map[Int, String], Map[String, Seq[Int]]) = {
      //iterate through each line and process only path or genome lines
      openFileWithIterator(gfa).foldLeft((Map[Int, String](), Map[String, Seq[Int]]())) {
        case ((graph, paths), line) => {
          if (line.startsWith("S")) {
            val (node, seq) = parseSegmentLineWithSeq(line)
            (graph + (node -> seq), paths)
          }
          else if (line.startsWith("P")) {
            //split line
            val split = line.split("\t")
            //assure genome has never been seen before
            assert(paths.get(split.head) == None)
            //add path to map
            (graph, paths + (split(1) -> split(2).split("""\W+""").toSeq.map(_.toInt)))
          } else (graph, paths)
        }
      }
    }
*/
    /**
      * Function to parse segment line
      * @return
      */
    def parseSegmentLine: String => Int = line => line.split("\t")(1).toInt

    def parseSegmentLineWithSeq: String => (Int,String) = line => {
      val split = line.split("\t")
      (split(1).toInt, split(2))
    }

    /**
      * Function to parse link line
      * @return
      */
    def parseLinkLine: String => (Int, Int) = line => {
      val split = line.split("\t")
      (split(1).toInt, split(3).toInt)
    }

  }

}
