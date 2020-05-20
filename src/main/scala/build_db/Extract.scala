package build_db

import java.io.{File, PrintWriter}

import atk.FastaIterator
import utilities.FileHandling.{tLines, timeStamp, verifyDirectory, verifyFile}
import utilities.{GFFutils, MinimapUtils}
import utilities.ConfigHandling

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
                     repeatRadius: Int = 10,
                     showWarnings: Boolean = false,
                     splitOverlaps: Double = 0.15,
                     minInterSize: Int = 11,
                     isCircular: Boolean = false,
                     useGene: Boolean = true,
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
      opt[Int]('i', "min-intergenic") action { (x, c) =>
        c.copy(minInterSize = x)
      } text ("Minimum size of an intergenic sequence to output to database (default is 11).")
      opt[Double]("split-overlaps") action { (x, c) =>
        c.copy(splitOverlaps = x)
      } text ("Split ORFs that overlap by this percentage using the smallest size of the two (default is 0.15).")
      opt[Unit]("circular") action { (x, c) =>
        c.copy(isCircular = true)
      } text ("Genomes are circular (turned off by default).")
      opt[Int]('r', "repeat-radius") action { (x, c) =>
        c.copy(repeatRadius = x)
      } text ("Remove any edge that are not within r positions away (default is 10). Used for identifying repetative" +
        " regions")
      opt[Unit]("use-cds") action { (x, c) =>
        c.copy(useGene = false)
      } text("Use 'CDS' annotations instead of 'gene' annotations (default is false, use 'gene' annotations).")
      opt[Unit]("show-warnings") action { (x, c) =>
        c.copy(showWarnings = true)
      }
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      }
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      println(config.outputDir)
      verifyDirectory(config.outputDir)
      verifyFile(config.genomesFile)
      extract(config)
    }
  }

  def extract(config: Config): Unit = {

	// DATA HANDLING

	// current configuration parameters
  	val params = ConfigHandling.getConfigParams(config)


	// call the function to store the parameters
	ConfigHandling.saveConfigParams("./testConfigWrite",params)

	// call the function to store the parameters
	val readParams = ConfigHandling.readConfigParams("./testConfigWrite")

	
	// get the new combined configuration file
	val config2 = ConfigHandling.generateConfig(readParams)

	println(config2)
	println(config)
	//

    //get annotation type
    val annotation_type = if(config.useGene) "gene" else "CDS"
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
    val pw_orf_id = new PrintWriter(config.outputDir + "/orf2id_mapping.txt")
    val pw_inter_id = new PrintWriter(config.outputDir + "/inter2id_mapping.txt")
    //output file for repetative regions
    val pw_repetative_regions = new PrintWriter(config.outputDir + "/repetative_regions.txt")
    //output file for id to fasta
    val pw_id2fasta = new PrintWriter(config.outputDir + "/id2fasta.txt")
    //curried gff parser based on parameters
    val gff_parser = parseORFs(config.isCircular, config.splitOverlaps, config.minInterSize, config.showWarnings) _

    /**
      * NOTE: created mutable variables since there is no need for parallalization and thus immutability. also to
      * improve legibility of code contrast to previous versions.
      **/
    // global orf id counter
    var global_orf_id = 0
    // global intergenic id counter
    var global_inter_id = 0
    // sample name -> chrm sequences name
    var map_y = scala.collection.mutable.HashMap[String, Seq[String]]()
    // sequence ID -> sorted orf IDS. Note: sequence ID is in format: <sample name>_<chrm sequence name>
    var map_z = scala.collection.mutable.HashMap[String, Seq[Int]]()
    var map_z_prime = scala.collection.mutable.HashMap[Int, String]()

    /**
      * Function to upate map z and z prime given the sequence ID and the ORF ID
      *
      * @return Unit
      */
    def updateZ: (String, Int) => Unit = (sequence_id, orf_id) => {
      //get current orfs
      val current_z = map_z.getOrElse(sequence_id, Seq[Int]())
      //update
      map_z += (sequence_id) -> (current_z.:+(orf_id))
      //sanity check
      assert(map_z_prime.get(orf_id).isEmpty, "Multiple occurrance of ORF " + orf_id)
      //update map
      map_z_prime += (orf_id -> sequence_id)
    }

    /**
      * Function to update map y and y prime given sample name and chrm sequence name
      *
      * @return Unit
      */
    def updateY: (String, String) => Unit = (sample_name, sequence_id) => {
      //update chrm sequences for current genome
      val current_y = map_y.getOrElse(sample_name, Seq[String]())
      //update map y
      map_y += (sample_name -> (current_y.:+(sequence_id)))
    }

    //iterate through each sample and construct database
    genomes.foreach { case (sample_name, assembly_file, gff_file) => {
      println(timeStamp + "--" + sample_name)
      //open GFF file, attempt to convert gff line to gff objects, collect non-empty, return those with proper
      // annotation type)
      val annotations = tLines(gff_file).map(x => toGFFLine(x, config.showWarnings))
        .collect{case x if(x.nonEmpty) => x.get}.filter(_.feature == annotation_type)
      if (config.verbose) println(timeStamp + "----Found " + annotations.size + " ORFs")
      //create output file to store ORF sequences
      val orf_sequences_file = new File(config.outputDir + "/" + sample_name + ".orfs.sequences.fasta")
      val pw_orf_seq = new PrintWriter(orf_sequences_file)
      val pw_inter_seq = new PrintWriter(config.outputDir + "/" + sample_name + ".intergenic.sequences.fasta")
      //open fasta assembly
      val assembly = new FastaIterator(assembly_file)
      //iterate through each chrm in the genome and output information to database
      while (assembly.hasNext) {
        //load fasta entry
        val current_chrm = assembly.next
        //size of contig
        val chrm_size = current_chrm.getSequence.size
        //get chrm name
        val chrm_name = current_chrm.getDescription.substring(1).split("\\s+").head
        val sequence_id = sample_name + "_" + chrm_name
        if(config.verbose) println(timeStamp + "----Processing sequence " + chrm_name)
        //get all (valid) corresponding annotations
        val corresponding_annotations =
          annotations.filter(x => x.chrm == chrm_name && x.start >= 1 && x.end <= chrm_size)
        //fetch all orfs and intergenic sequences
        val (orfs, updated_orf_id, intergenics, updated_inter_id) = gff_parser(current_chrm, sample_name,
          global_orf_id, global_inter_id, corresponding_annotations)
        if(config.verbose) println(timeStamp + "----Curated " + orfs.size + " ORFs and " + intergenics.size +
          " intergenic sequences")
        //output orf2id mapping scheme and orf sequence
        orfs.foreach(orf => {
          //update map z and z prime
          updateZ(sequence_id, orf.id)
          //output to database
          pw_id2fasta.println(orf.id + "\t" + sample_name + ".orfs.sequences.fasta")
          pw_orf_id.println(orf.description + "\t" + sample_name + "\t" + orf.id)
          pw_orf_seq.println(">" + orf.id + "\n" + current_chrm.getSequence.substring(orf.start - 1, orf.end))
        })
        //update global orf id
        global_orf_id = updated_orf_id
        //output inter2id mapping scheme and intergenic sequence
        intergenics.foreach(inter => {
          //update global intergenic id
          global_inter_id = inter.id
          pw_inter_id.println(inter.description + "\t" + inter.id)
          //note that there can be multiple linear intergenic coordinates if the genome is circular
          pw_inter_seq.println(">" + inter.id + "\n" +
            inter.coords.map(x => current_chrm.getSequence.substring(x._1 - 1, x._2)).mkString(""))
        })
        //update map y and y prime
        updateY(sample_name, sequence_id)
        //update global intergenic id
        global_inter_id = updated_inter_id
      }
      //close output files
      pw_orf_seq.close()
      pw_inter_seq.close()

      //get repetative region start,end node positions
      val repetative_regions = {
        //run minimap with default parameters and no self alignments turned on
        val rgraph = {
          collectBestAlignments(orf_sequences_file, orf_sequences_file, 0.75, 11, 5, HashMap[Int, Set[Int]](), true)
            //map,reduce approach to obtain repetative regions
            .map { case (node, edges) => {
            //gets seq ID for node
            val node_seq = map_z_prime(node)

            //function to compute distance of two nodes
            def dist: Int => Int = e => Math.abs(e - node)
            //remove any edge whose node is not in the same sequence and within specified repeat radius
            (node, edges.filter(edge => map_z_prime(edge) == node_seq && dist(edge) > 0 &&
              dist(edge) <= config.repeatRadius))
          }
          }.filter(!_._2.isEmpty)
            //construct undirected version of graphs
            .foldLeft(HashMap[Int, Set[Int]]()) { case (g, (node, edges)) => {
            edges.foldLeft(g)((local_g, edge) => {
              val current_edge = local_g.getOrElse(edge, Set[Int]())
              local_g + (edge -> (current_edge + node))
            }) + (node -> edges)
          }
          }
        }

        if (config.verbose) println(timeStamp + "--Found repeat graph of at least " + rgraph.size + " nodes")

        /**
          * Method to find connected components in the repetative regions graph using breadth-first search traversal
          *
          * @param nodes
          * @param ccs
          * @return
          */
        def getRepeatRegionsCC(nodes: Set[Int], ccs: List[Set[Int]]): List[Set[Int]] = {
          /**
            * BFS traversal
            *
            * @param toVisit      Nodes to visit, in-order
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

          if (nodes.isEmpty) ccs
          else {
            //get connected component for current node
            val cc = breadthFirstSearch(Seq(nodes.head), Set())
            //remove all nodes in the connectec component and move on
            getRepeatRegionsCC(nodes.tail -- cc, cc :: ccs)
          }
        }
        //get repeat regions CCs in the repeat graph, then summarize each CC as (min,max) node
        getRepeatRegionsCC(rgraph.flatMap(x => x._2 + (x._1)).toSet, List()).map(_.toSeq.sortBy(identity))
      }

      /**
        * Function to determine whether two repeat region graphs overlap
        *
        * @return Boolean
        */
      def rangeOverlap: ((Int, Int), (Int, Int)) => Boolean = (r, q) => r._1 <= q._2 && r._2 >= q._1

      /**
        * Tail-recursive method to merge repeat regions that overlap in terms of node location
        *
        * @param regions
        * @param merged
        * @return
        */
      def mergeRepeatRegions(regions: List[(Int, Int)], merged: List[(Int, Int)]): List[(Int, Int)] = {
        regions match {
          case Nil => merged
          case head :: tail => {
            val overlaps = head :: tail.filter(g => rangeOverlap(head, g))
            val maximal_region = (overlaps.map(_._1).min, overlaps.map(_._2).max)
            mergeRepeatRegions(tail.filterNot(x => overlaps.contains(x)), maximal_region :: merged)
          }
        }
      }
      //output to file
      repetative_regions.foreach(x => pw_repetative_regions.println(x.mkString(",")))
    }}

    //output map y and y prime to database
    map_y.foreach{case(sample_name, sequence_ids) => {
      pw_Y.println(sample_name + "\t" + sequence_ids.mkString(","))
      sequence_ids.foreach(contig => pw_Y_prime.println(contig + "\t" + sample_name))
    }}

    //output map z and z prime to database
    map_z.foreach{case (sequence_id, orf_ids) => {
      pw_Z.println(sequence_id + "\t" + orf_ids.mkString(","))
      orf_ids.foreach(orf => pw_Z_prime.println(orf + "\t" + sequence_id))
    }}

    //close output files
    List(pw_Y, pw_Y_prime, pw_Z, pw_Z_prime, pw_orf_id, pw_repetative_regions, pw_id2fasta,
      pw_inter_id).foreach(_.close())
    println(timeStamp + "Successfully completed!")
  }


  /**
    * Function to determine features of an ORF. If it's a single orf, use as is. If not, there are overlapping ORFs
    * so artificially create an ORF with the max boundaries (smallest start, largest end)
    *
    * @return
    */
  def getMaximalORF: (List[GFFLine], Boolean) => (String, Int, Int) = (gffs, warning) => {
    if (gffs.size > 1) {
      //get left and right-most coordinates
      val min_start = gffs.map(_.start).min
      val max_end = gffs.map(_.end).max
      //log message
      if (warning) println(timeStamp + "----WARNING: overlapping ORFs detected. Using the following maximal start,end positions: " +
        (min_start, max_end))
      (makeGenericORFname(Option((max_end - min_start) + 1))(gffs.minBy(_.start)), min_start, max_end)
    } else (makeGenericORFname(None)(gffs.head), gffs.head.start, gffs.head.end)
  }

  /**
    * Function to create a generic name for a GFF line
    *
    * @return String concatenation of (chrm,name,start position)
    */
  def makeGenericORFname(size: Option[Int]): GFFLine => String = gff => {
    if (gff.name == None) Seq(gff.chrm, gff.start).mkString("_")
    else Seq(gff.chrm, gff.name.get, gff.start, if (!size.isEmpty) size.get else (gff.end - gff.start) + 1).mkString("_")
  }

  /**
    * Function to parse genomes file
    *
    * @return Tuple of (Genome ID, FASTA assembly path, GFF file path)
    */
  def getGenomesInfo: String => (String, File, File) = line => {
    val tmp = line.split("\t")
    assume(tmp.size >= 3, "Unexpected number of entries provided for the following line: " + line)
    (tmp.head, new File(tmp(1)), new File(tmp(2)))
  }

  /**
    * Function to verify that the assemblies and gff files are valid paths
    *
    * @return None (unit)
    */
  def verifyCorrectGenomeEntries: (String, File, File) => Unit = (genome_id, assembly, gff) => {
    assume(assembly.exists() && assembly.isFile, "Invalid assembly file for " + genome_id + ": " + assembly
      .getAbsolutePath)
    assume(gff.exists() && gff.isFile, "Invalid GFF file for " + genome_id + ": " + assembly.getAbsolutePath)
  }

  /**
    * Function to verify that a collection of any type is unique
    *
    * @return Boolean
    */
  def isUnique: Iterable[Any] => Unit = collection =>
    assume(collection.toList.distinct.size == collection.toList.size, "The genome ID's provided are not all unique")

}
