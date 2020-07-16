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
    /**
      * Load file containing a map from ORF id to Node Id
      * @return Map[Int,Int]
      */
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
      _orfid2nodeid
    }

    /**
      * Load file containing a map from Node ID to list of gene names
      * @return Map[Int,List[String]]
      */
    def loadNodeID2gene: File => scala.collection.mutable.Map[Int,List[String]] = db => {
      // get file from the database that contains the desired mapping
      val file = db.listFiles().find(_.getName == "orf2node_id.txt")
      // check if file really exists before parsing
      assert(file != None, "Could not find file mapping orf ID to node ID.")
      // iterate through file and create map
      val _nodeID2gene = openFileWithIterator(file.get).
          foldLeft(scala.collection.mutable.Map[Int,List[String]]())((map, line) => {
            // delimiter are tabs
            val split = line.split("\t")
            // get gene name 
            val gene = split(0).toString
            // get Node Identifiers
            val nodeID = split(3).toInt
            // append to the map
            map.get(nodeID) match {
              case Some(xs:List[String]) => map.updated(nodeID, xs :+ gene)
              case None => map += (nodeID -> List(gene))
            }
          })
      _nodeID2gene
    }

    /**
      * Load file containing a mapping from the Node ID to list of strains
      * @return Map[Int,List[String]]
      */
    def loadNodeID2strain: File => scala.collection.mutable.Map[Int,List[String]] = db => {
      // get file from the database that contains the desired mapping
      val file = db.listFiles().find(_.getName == "orf2node_id.txt")
      // check if file really exists before parsing
      assert(file != None, "Could not find file mapping orf ID to node ID.")
      // iterate through file and create map
      val _nodeID2strain = openFileWithIterator(file.get).
          foldLeft(scala.collection.mutable.Map[Int,List[String]]())((map, line) => {
            // delimiter are tabs
            val split = line.split("\t")
            // get strain name 
            val strain = split(1).toString
            // get Node Identifiers
            val nodeID = split(3).toInt
            // append to the map
            map.get(nodeID) match {
              case Some(xs:List[String]) => map.updated(nodeID, xs :+ strain)
              case None => map += (nodeID -> List(strain))
            }
          })
      _nodeID2strain.foreach{
                case (key, value) => _nodeID2strain(key) = value.sorted}
      _nodeID2strain
    }

    /**
      * From database, load the unique strains name
      * @return List(String)
      */
    def loadStrainNames: File => List[String] = db => {
      // get file from those of the database that contains the desired mapping
      val file = db.listFiles().find(_.getName == "orf2node_id.txt")
      // check if file really exists before parsing
      assert(file != None, "Could not find file mapping orf ID to node ID.")
      // iterate through file and create a seq with strain names
      val _strains = openFileWithIterator(file.get).
          foldLeft(List[String]())((seq, line) => {
          // get strain name field 
          val strain = line.split("\t")(1).toString
          // append to the list (even if it exists)
          seq :+ strain 
        })
      // delete duplicates and sort names
      _strains.distinct.sorted
    }

    /**
      * From database, load the unique strains and generate strain intersections 
      * @return Map[List(String),Int]
      */
    def strainIntersection: File => Map[List[String],Int] = db => {
      // get file from those of the database that contains the desired mapping
      val file = db.listFiles().find(_.getName == "orf2node_id.txt")
      // check if file really exists before parsing
      assert(file != None, "Could not find file mapping orf ID to node ID.")
      // iterate through file and create a seq with strain names
      val _strains = openFileWithIterator(file.get).
          foldLeft(scala.collection.mutable.Seq[String]())((seq, line) => {
          // get strain name field 
          val strain = line.split("\t")(1).toString
          // append to the list (even if it exists)
          seq :+ strain 
        })
      // delete duplicates, generate all possible subsets and convert to map
      val output = _strains.toSet[String].subsets.map(subset => subset.toList.sorted -> 0).toMap
      // remove the empty set and return
      output - List()
    }
  }

  /**
    * Tools for using Ptolemy's graph index
    */
  trait GraphIndex {
    /**
      * Load global node kmer index in the format of Map(kmer hashvalue, Set(node IDs))
      *
      * @return Map[Int, Set[Int]
      **/
    def loadGlobalNodeKmerIndex: File => Map[Int, Set[Int]] = file => {
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
      * Load global intergenic kmer index in the format of Set[hashvalue]
      * @return Set[Int]
      */
    def loadGlobalInterKmerIndex: File => Set[Int] = file => {
      //iterate through eac line in the file
      openFileWithIterator(file).foldLeft(Set[Int]())((kmap, line) => kmap + (line.toInt))
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
