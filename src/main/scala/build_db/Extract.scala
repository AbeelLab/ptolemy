package build_db

import java.io.{File, PrintWriter}

import atk.FastaIterator
import utilities.FileHandling.{tLines, timeStamp, verifyDirectory, verifyFile}
import utilities.GFFutils
import utilities.MinimapUtils

import scala.annotation.tailrec
import scala.collection.immutable.HashMap

/**
  * Author: Alex N. Salazar
  * Created on 14-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object Extract extends tLines with GFFutils with MinimapUtils {

  case class Config(
                     genomesFile: File = null,
                     outputDir: File = null,
                     singleHeader: Boolean = false,
                     repeatRadius: Int = 10,
                     showWarnings: Boolean = false,
                     verbose: Boolean = false
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("extract") {
      opt[File]('g', "genomes") required() action { (x, c) =>
        c.copy(genomesFile = x)
      } text ("Tab-delimited file containing (genome ID, FASTA path, and GFF file), on per line.")
      opt[File]('o', "output-database") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory for database to be stored.")
      note("\nOPTIONAL\n")
      opt[Int]('r', "repeat-radius") action { (x, c) =>
        c.copy(repeatRadius = x)
      } text ("Remove any edge that are not within r positions away (default is 10). Used for identifying repetative" +
        " regions")
      opt[Unit]("single-header") action { (x, c) =>
        c.copy(singleHeader = true)
      } text ("Split FASTA sequence names by whitespace and use the first element as the sequence name rather than " +
        "using the whole string (turned off by default).")
      opt[Unit]("show-warnings") action { (x, c) =>
        c.copy(showWarnings = true)
      }
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      }
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.genomesFile)
      extract(config)
    }
  }

  def extract(config: Config): Unit ={
    //open list of genomes
    val genomes = tLines(config.genomesFile).map(getGenomesInfo(_))
    println(timeStamp + "Found " + genomes.size + " genome entries")
    println(timeStamp + "--Verifying paths")
    //verify that all paths provided are valid
    genomes.foreach(x => verifyCorrectGenomeEntries(x._1, x._2, x._3))
    //verify that all genome IDs are unique
    isUnique(genomes.map(_._1))
    //create global output files for storing hashmap
    val pw_Y = new PrintWriter(config.outputDir + "/global_y.txt")
    val pw_Y_prime = new PrintWriter(config.outputDir + "/global_y_prime.txt")
    val pw_Z = new PrintWriter(config.outputDir + "/global_z.txt")
    val pw_Z_prime = new PrintWriter(config.outputDir + "/global_z_prime.txt")
    //create output file to store mapping between ORF and assigned ORF ID
    val pw_id_mapping = new PrintWriter(config.outputDir + "/orf2id_mapping.txt")
    //output file for repetative regions
    val pw_repetative_regions = new PrintWriter(config.outputDir + "/repetative_regions.txt")
    //output file for id to fasta
    val pw_id2fasta = new PrintWriter(config.outputDir + "/id2fasta.txt")
    /**
      * Method to construct database for each given genome
      * @param genome_id Starting integer to be used as global ORF ID (will be incremented
      * @param fasta File object of FASTA assembly
      * @param gff File object of GFF file
      * @param dir
      * @param global_orf_id
      * @return
      */
    def constructDataBase(genome_id: String, fasta: File, gff: File, dir: File, global_orf_id: Int): Int = {

      //open GFF file
      val orfs = tLines(gff).map(toGFFLine(_)).filter(_.feature == "gene")
      if(config.verbose) println(timeStamp + "----Found " + orfs.size + " ORFs")
      val sequences_file = new File(dir + "/" + genome_id + ".orfs.sequences.fasta")
      //create output file to store ORF sequences
      val pw_sequences = new PrintWriter(sequences_file)
      //open fasta assembly
      val assembly = new FastaIterator(fasta)

      /**
        * Internal function to output information to database given a set of ORFs, the global ID, and the fasta entry
        * @return Unit
        */
      def outputToDB: (List[GFFLine], Int, atk.Record, Int) => Unit = (overlaps, id, fasta_entry, size) => {
        //get maximal ORF information
        val (generic_name, start, end) = getMaximalORF(overlaps, config.showWarnings)
        if(end > size) {
          if(config.showWarnings) println(timeStamp + "--WARNING: Skipping ORF with coordinates outside of sequence " +
            "boundary: " + (generic_name, start, end))
        }
        else {
          //output generic name and unique ID assigned
          pw_id_mapping.println(Seq(generic_name, genome_id, id).mkString("\t"))
          pw_id2fasta.println(Seq(id, genome_id + ".orfs.sequences.fasta").mkString("\t"))
          //get sequence
          val sequence = fasta_entry.getSequence.substring(start - 1, end)
          //output sequence
          pw_sequences.println(Seq(">" + id, sequence).mkString("\n"))
          if (config.verbose) print(id + ",")
        }
      }

      /**
        * Method to iterate through each FASTA assembly and extract:
        *   -All sequences (contigs/chrms) in assembly (Y hashmap)
        *   -All ORFs in a sequence (Z hashmap)
        *   -Assign each ORF an unique global ID
        *   -Output ORF mapping and sequences to database
        * @param orf_id Integer to to start with for global ID
        * @param Y Hashmap of genome -> sequences
        * @param Z Hashmap of sequence -> ORFS (ordered)
        * @return (New starting ID as INT, Y, Z)
        */
      def parseORFs(orf_id: Int,
                    Y: HashMap[String, Seq[String]],
                    Z: HashMap[String, Seq[Int]]): (Int, HashMap[String, Seq[String]], HashMap[String, Seq[Int]]) = {
        if (!assembly.hasNext) (orf_id, Y, Z)
        else {
          //get current fasta entry
          val current_entry = assembly.next()
          val entry_size = current_entry.getSequence.size
          //get id
          val entry_id = {
            if(config.singleHeader) current_entry.getDescription.substring(1).split("\\s+").head
            else current_entry.getDescription.substring(1)
          }
          val entry_id_and_genome_id = genome_id + "_" + entry_id
          //find corresponding ORFs in sequence
          val corresponding_orfs = orfs.filter(x => x.chrm == entry_id)
          if(config.verbose) println(timeStamp + "------Found " + corresponding_orfs.size + " ORFs in " + entry_id +
            "\n------Adding:")
          //iterate through each orf and output to database
          val (updated_Z, last_orf_ID, last_overlaps) = corresponding_orfs.foldLeft((Z, orf_id, List[GFFLine]())){
            case((local_Z,local_orf_id,overlaps), current_orf) => {
              //handle case of first iteration and cases of overlapping ORFs
              if (overlaps.isEmpty || overlaps.exists(isOverlap(_, current_orf)))
                (local_Z, local_orf_id, current_orf :: overlaps)
              //no overlapping annotations
              else {
                //output information to database
                outputToDB(overlaps, local_orf_id, current_entry, entry_size)
                //get current sequence of ORFs
                val current = local_Z.getOrElse(entry_id_and_genome_id, Seq[Int]())
                //update local hashmap Z and the global unique ORF identifier
                (local_Z + (entry_id_and_genome_id -> current.:+(local_orf_id)),
                  local_orf_id + 1,
                  List(current_orf))
              }
            }
          }
          //updated Y hashmap
          val updated_Y = {
            //get current sequences in genome
            val current = Y.getOrElse(genome_id, Seq[String]())
            //update sequence to genome
            Y + (genome_id -> current.:+(entry_id_and_genome_id))
          }
          //no overlapping orfs at the end of a sequence
          if(last_overlaps.isEmpty) {
            if(config.verbose)println
            parseORFs(last_orf_ID, updated_Y, updated_Z)
          }
          //one last set of overlapping orfs at end of sequence
          else {
            //output to DB
            outputToDB(last_overlaps, last_orf_ID, current_entry, entry_size)
            //get current sequence of ORFs
            val current = updated_Z.getOrElse(entry_id_and_genome_id, Seq[Int]())
            if(config.verbose)println
            //update
            parseORFs(last_orf_ID+1, updated_Y, updated_Z + (entry_id_and_genome_id -> current.:+(last_orf_ID)))
          }
        }
      }
      val (new_id, y_hashmap, z_hashmap) = parseORFs(global_orf_id, HashMap(), HashMap())
      //create hashmap of orf -> sequence id
      val z_hashmap_prime = z_hashmap.foldLeft(HashMap[Int,String]())((h,s) =>
        s._2.foldLeft(h)((local_h,o) => local_h + (o -> s._1)))
      //update global Y and Y'
      y_hashmap.foreach{case(genome, sequences) => {
        pw_Y.println(Seq(genome, sequences.mkString(",")).mkString("\t"))
        sequences.foreach(sequence => pw_Y_prime.println(Seq(sequence, genome).mkString("\t")))
        sequences.foreach(sequence => pw_Y_prime.println(Seq(sequence, genome).mkString("\t")))
      }}
      //update global Z and Z'
      z_hashmap.foreach{case(sequence, orfs) => {
        pw_Z.println(Seq(sequence, orfs.mkString(",")).mkString("\t"))
        orfs.foreach(orf => pw_Z_prime.println(Seq(orf, sequence).mkString("\t")))
      }}
      pw_sequences.close
      //get repetative region start,end node positions
      val repetative_regions = {
        //run minimap with default parameters and no self alignments turned on
        val rgraph = collectBestAlignments(sequences_file, sequences_file, 0.75, 11, 5, HashMap[Int, Set[Int]](), true)
          //map,reduce approach to obtain repetative regions
              .map{case(node,edges) =>{
                //gets seq ID for node
                val node_seq = z_hashmap_prime(node)
                val dist: Int => Int = e => Math.abs(e - node)
                //remove any edge whose node is not in the same sequence and within specified repeat radius
                (node, edges.filter(edge => z_hashmap_prime(edge) == node_seq && dist(edge) > 0 &&
                  dist(edge) <= config.repeatRadius))
              }}.filter(!_._2.isEmpty)
                //construct undirected version of graphs
                .foldLeft(HashMap[Int,Set[Int]]()){case (g, (node,edges)) => {
                  edges.foldLeft(g)((local_g, edge) => {
                    val current_edge = local_g.getOrElse(edge, Set[Int]())
                    local_g + (edge -> (current_edge + node))
                  }) + (node -> edges)
                }}
        if(config.verbose) println(timeStamp + "--Found repeat graph of at least " + rgraph.size + " nodes")
        /**
          * Method to find connected components in the repetative regions graph using breadth-first search traversal
          * @param nodes
          * @param ccs
          * @return
          */
        def getRepeatRegionsCC(nodes: Set[Int], ccs: List[Set[Int]]): List[Set[Int]] = {
          /**
            * BFS traversal
            * @param toVisit Nodes to visit, in-order
            * @param localVisited Nodes visited in the current traversal
            * @return
            */
          @tailrec def breadthFirstSearch(toVisit: Seq[Int], localVisited: Set[Int]): Set[Int] = {
            toVisit match {
              case Nil => localVisited
              case head :: tail => {
                //get neighbours of current node
                val neighbours = rgraph.getOrElse(head, Set[Int]()).toSeq.filterNot(localVisited.contains(_))
                //traverse neighbours
                breadthFirstSearch(tail ++ neighbours, localVisited + head)
              }
            }
          }
          if(nodes.isEmpty) ccs
          else {
            //get connected component for current node
            val cc = breadthFirstSearch(Seq(nodes.head), Set())
            //remove all nodes in the connectec component and move on
            getRepeatRegionsCC(nodes.tail -- cc, cc :: ccs)
          }
        }
        //get repeat regions CCs in the repeat graph, then summarize each CC as (min,max) node
        getRepeatRegionsCC(rgraph.flatMap(x => x._2 +(x._1)).toSet, List()).map(_.toSeq.sortBy(identity))
      }

      /**
        * Function to determine whether two repeat region graphs overlap
        * @return Boolean
        */
      def rangeOverlap: ((Int,Int), (Int,Int)) => Boolean = (r,q) => r._1 <= q._2 && r._2 >= q._1

      /**
        * Tail-recursive method to merge repeat regions that overlap in terms of node location
        * @param regions
        * @param merged
        * @return
        */
      def mergeRepeatRegions(regions: List[(Int,Int)], merged: List[(Int,Int)]): List[(Int,Int)] = {
        regions match {
          case Nil => merged
          case head :: tail => {
            val overlaps = head :: tail.filter(g => rangeOverlap(head,g))
            val maximal_region = (overlaps.map(_._1).min, overlaps.map(_._2).max)
            mergeRepeatRegions(tail.filterNot(x => overlaps.contains(x)), maximal_region :: merged)
          }
        }
      }
      //val maximal_repetative_regions = mergeRepeatRegions(repetative_regions, List())
      //output to file
      repetative_regions.foreach(x => pw_repetative_regions.println(x.mkString(",")))
      new_id
    }
    println(timeStamp + "Constructing database: ")
    //iterate through each assembly and construct database
    genomes.foldLeft(0){ case (id,genome) => {
      println(timeStamp + "--" + genome._1)
      constructDataBase(genome._1, genome._2, genome._3, config.outputDir, id)
    }}
    //close output files
    List(pw_Y, pw_Y_prime, pw_Z, pw_Z_prime, pw_id_mapping, pw_repetative_regions, pw_id2fasta).foreach(_.close())
    println(timeStamp + "Successfully completed!")
  }



  def getMaximalORF: (List[GFFLine], Boolean) => (String,Int,Int) = (gffs, warning) => {
    val largest = gffs.maxBy(x => (x.end - x.start)+1)
    if(gffs.size > 1 && warning) {
      println(timeStamp + "----WARNING: overlapping ORFs detected. Using the following maximal start,end positions: " +
        (largest.start,largest.end))
    }
    (makeGenericORFname(largest), largest.start, largest.end)
  }

  /**
    * Function to create a generic name for a GFF line
    * @return String concatenation of (chrm,name,start position)
    */
  def makeGenericORFname: GFFLine => String = gff => {
   if(gff.name == None) Seq(gff.chrm, gff.start).mkString("_")
   else Seq(gff.chrm, gff.name.get, gff.start, (gff.end-gff.start)+1).mkString("_")
  }

  /**Method to determine whether two GFF objects overlap*/
  def isOverlap(x: GFFLine, y: GFFLine): Boolean = x.start <= y.end && x.end >= y.start

  /**
    * Function to parse genomes file
    * @return Tuple of (Genome ID, FASTA assembly path, GFF file path)
    */
  def getGenomesInfo: String => (String,File,File) = line => {
    val tmp = line.split("\t")
    assume(tmp.size >= 3, "Unexpected number of entries provided for the following line: " + line)
    (tmp.head, new File(tmp(1)), new File(tmp(2)))
  }

  /**
    * Function to verify that the assemblies and gff files are valid paths
    * @return None (unit)
    */
  def verifyCorrectGenomeEntries: (String, File, File) => Unit = (genome_id, assembly, gff) => {
    assume(assembly.exists() && assembly.isFile, "Invalid assembly file for " + genome_id + ": " + assembly
      .getAbsolutePath)
    assume(gff.exists() && gff.isFile, "Invalid GFF file for " + genome_id + ": " + assembly.getAbsolutePath)
  }

  /**
    * Function to verify that a collection of any type is unique
    * @return Boolean
    */
  def isUnique: Iterable[Any] => Unit = collection =>
    assume(collection.toList.distinct.size == collection.toList.size, "The genome ID's provided are not all unique")

}
