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
      * Load canonical quiver from a GFA-formatted file. Constructs the canonical quiver graph, all genome paths in
      * th graph, all genomes in the graph, and the node coverages
      *
      * @param gfa GFA file
      * @return Graph as Map[Int, Set[Int], Paths as Map[String, Seq[Int], genomes as Map[String, Set[String], and
      *         node coverages as Map[Int, Int]
      */
    def loadGFA(gfa: File): (Map[Int, Set[Int]], Map[String, Seq[Int]], Map[String, Set[String]], Map[Int, Int]) = {
      //iterate through each line and construct graph, paths, and genome data structures
      val (graph, paths, genomes) = openFileWithIterator(gfa)
        .foldLeft((Map[Int, Set[Int]](), Map[String, Seq[Int]](), Map[String, Set[String]]())) {
          case ((graph, paths, genomes), line) => {
            //node line
            if (line.startsWith("S")) {
              //parse node
              val node = parseSegmentLine(line)
              //sanity check: node line should only be seen once
              assert(graph.get(node).isEmpty, "Observed segment line for a node multiple times: " + line)
              //update
              (graph + (node -> Set[Int]()), paths, genomes)
            }
            //edge line
            else if (line.startsWith("L")) {
              //parse node and corresponding edge
              val (node, edge) = parseLinkLine(line)
              //get current edges
              val current = graph.getOrElse(node, Set[Int]())
              //update
              (graph + (node -> (current + edge)), paths, genomes)
            }
            //genome line
            else if (line.startsWith("G")) {
              //split line
              val split = line.split("\t")
              //sanity check: genome line for current genome should be observed only once
              assert(genomes.get(split(1)) == None, "Multiple occurrances of a genome line for a genome: " + line)
              //update
              (graph, paths, genomes + (split(1) -> split(2).split(",").toSet))
            }
            //path line
            else if (line.startsWith("P")) {
              //split line
              val split = line.split("\t")
              //assure genome has never been seen before
              assert(paths.get(split.head) == None, "Multiple occurrances of path line for a genome: " + line)
              //add path to map
              (graph, paths + (split(1) -> split(2).split("""\W+""").toSeq.map(_.toInt)), genomes)
            }
            //some other line, skip
            else (graph, paths, genomes)
          }
        }
      //obtain coverages of every node in the graph
      val coverages = paths.values.foldLeft(Map[Int, Int]())((node_coverages, path) => {
        //iterate through each path and update coverage of node
        path.foldLeft(node_coverages)((local_coverage, node) => {
          //get current coverage
          val current = local_coverage.getOrElse(node, 0)
          //increment coverage
          local_coverage + (node -> (current + 1))
        })
      })
      //return
      (graph, paths, genomes, coverages)
    }

    /**
      * Function to load only the nodes in the GFA
      *
      * @param gfa GFA file
      * @return Set[Int]
      */
    def loadNodesGFA(gfa: File): Set[Int] = {
      //iterate through each line and process only segment lines
      openFileWithIterator(gfa).foldLeft(Set[Int]())((nodes, line) => {
        if (!line.startsWith("S")) nodes else nodes + (parseSegmentLine(line))
      })
    }

    /**
      * Function to parse path alignment line. Returns 2-tuple of (read ID, ORFs in order of the read)
      *
      * @return (String, IndexedSeq[Int])
      */
    def parseAlignment: String => (String, IndexedSeq[Int]) = line => {
      val split = line.split("\t")
      (split(1), if (split(2).isEmpty) IndexedSeq[Int]() else split(2).split("""\W+""").map(_.toInt).toIndexedSeq)
    }

    /**
      * Get only the graph and the paths in a GFA file.
      *
      * @return Map as ID -> Seq[Node IDs]
      *         *
      *         def getSequencesPaths(gfa: File): (Map[Int, String], Map[String, Seq[Int]) = {
      *         //iterate through each line and process only path or genome lines
      *         openFileWithIterator(gfa).foldLeft((Map[Int, String](), Map[String, Seq[Int]())) {
      *         case ((graph, paths), line) => {
      *         if (line.startsWith("S")) {
      *         val (node, seq) = parseSegmentLineWithSeq(line)
      *         (graph + (node -> seq), paths)
      *         }
      *         else if (line.startsWith("P")) {
      *         //split line
      *         val split = line.split("\t")
      *         //assure genome has never been seen before
      *         assert(paths.get(split.head) == None)
      *         //add path to map
      *         (graph, paths + (split(1) -> split(2).split("""\W+""").toSeq.map(_.toInt)))
      *         } else (graph, paths)
      *         }
      *         }
      *         }
      */
    /**
      * Function to parse segment line and return node
      *
      * @return Int
      */
    def parseSegmentLine: String => Int = line => line.split("\t")(1).toInt

    def parseSegmentLineWithSeq: String => (Int, String) = line => {
      val split = line.split("\t")
      (split(1).toInt, split(2))
    }

    /**
      * Function to parse link line and return edge as 2-tuple (node, node)
      *
      * @return (Int,Int)
      */
    def parseLinkLine: String => (Int, Int) = line => {
      val split = line.split("\t")
      (split(1).toInt, split(3).toInt)
    }

  }

}
