package utilities

import java.io.File
import utilities.MinimizerUtils.{Minimizer, MMethods}
import utilities.FileHandling.{openFileWithIterator}

/**
  * Author: Alex N. Salazar
  * Created on 22-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object DatabaseUtils extends MMethods {

  trait PtolemyDB {
    //create map from orf ID to node IDs
    def loadORFid2Nodeid: File => Map[Int,Int] = db => {
      //fetch file containing mapping of orf id to node id
      val file = db.listFiles().find(_.getName == "orf2node_id.txt")
      //sanity check
      assert(file != None, "Could not find file mapping orf ID to node ID.")
      //iterate through file and create map
      val _orfid2nodeid = openFileWithIterator(file.get).foldLeft(Map[Int, Int]())((map, line) => {
        val split = line.split("\t")
        //get IDs
        val (orf, node) = (split(2).toInt, split(3).toInt)
        //sanity check
        assert(map.get(orf) == None, "ORF ID " + orf + " appears multiple times at line " + line)
        map + (orf -> node)
      })
      //get all found node IDs
      val found_nodeids = _orfid2nodeid.map(_._2).toSet
      _orfid2nodeid
    }
  }

  /**
    * Tools for using Ptolemy's graph index
    */
  trait GraphIndex {
    /**
      * Load global kmer index in the format of Map(kmer hashvalue, Set(node IDs))
      *
      * @return Map[Int, Set[Int]
      **/
    def loadGlobalKmerIndex: File => Map[Int, Set[Int]] = file => {
      //iterate through each line in the file
      openFileWithIterator(file).foldLeft(Map[Int, Set[Int]]())((kmap, line) => {
        val line_split = line.split("\t")
        //first element is assumed to be the kmer hashvalue>
        val kmer = line_split.head.toInt
        //sanity check
        assert(kmap.get(kmer) == None, "Multiple occurances of kmer hash value: " + kmer + " line: " + line)
        //construct kmer's node set, every entry is a node id
        kmap + (kmer -> line_split(1).split(",").foldLeft(Set[Int]())((node_set, node_id) => node_set + (node_id.toInt)))
      })
    }

    /**
      * Load global node map in the format of Map[node id, MinKmer]
      *
      * @return
      */
    def loadGlobalNodeIndex: File => Map[Int, Seq[Minimizer]] = file => {
      // iterate through each line in the file
      openFileWithIterator(file).foldLeft(Map[Int, Seq[Minimizer]]())((nmap, line) => {
        //split line
        val tmp = line.split(":")
        //first element assumed to be node id
        val node_id = tmp.head.toInt
        //starting positions
        val minimizer = tmp(1).split(",").map(x => {
          val tmp = x.split("_")
          new Minimizer(tmp.head.toInt, tmp(1).toInt, tmp(2).toInt)
        }).toSeq
        assert(nmap.get(node_id) == None, "Multiple occurrances of node " + node_id + ": " + line)
        nmap + (node_id -> minimizer)
      })
    }

    /**
      * Construct global node info index in the format of Map[node id, (avg minimizer set size, avg seq size)]
      *
      * @return Map[Int, (Int, Int)]
      */
    def loadGlobanNodeInfoIndex: File => Map[Int, (Int, Int)] = file => {
      //iterate through each line in the file
      openFileWithIterator(file).foldLeft(Map[Int, (Int, Int)]())((nimap, line) => {
        //split line
        val tmp = line.split("\t")
        //get node id, average size of minimizer set, and average sequence size
        val (node_id, minimizer_set_size, seq_size) = (tmp(0).toInt, tmp(1).toInt, tmp(2).toInt)
        assert(nimap.get(node_id) == None, "Multiple occurrances of node " + node_id + ": " + line)
        nimap + (node_id -> (minimizer_set_size, seq_size))
      })
    }

  }

}
