package syntenic_anchors

import java.io.{File, PrintWriter}

import utilities.FileHandling._
import utilities.NumericalUtils.min
import utilities.{GFFutils, MinimapUtils}
import utilities.GraphUtils.getConnectedComponents
import atk.Tool.progress
import scala.collection.immutable.HashMap
import scala.collection.parallel.ForkJoinTaskSupport


/**
  * Author: Alex N. Salazar
  * Created on 14-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SyntenicAnchorsC extends GFFutils with MinimapUtils {

  case class Config(
                     database: File = null,
                     outputDir: File = null,
                     gamma: Double = 0.75,
                     flankingWindow: Int = 10,
                     verbose: Boolean = false,
                     mergeCC: Boolean = false,
                     dump: Boolean = false,
                     syntenicFraction: Double = 0.5,
                     kmerSize: Int = 11,
                     minimizerWindow: Int = 3,
                     isCircular: Boolean = false,
                     mergeSingleBRH: Boolean = false,
                     maxThreads: Int = 1,
                     brh: File = null
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("syntenic-anchors") {
      opt[File]("db") required() action { (x, c) =>
        c.copy(database = x)
      } text ("Directory path of database (e.g. output directory of 'extract' module).")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory.")
      note("\nOPTIONAL: SYNTENY PARAMETERS\n")
      opt[Unit]("circular") action { (x, c) =>
        c.copy(isCircular = true)
      } text ("Genome can be circular (default is false).")
      opt[Int]('f', "flanking-window") action { (x, c) =>
        c.copy(flankingWindow = x)
      } text ("Flanking window (default is 10).")
      opt[Unit]("merge-ccs") action { (x, c) =>
        c.copy(mergeCC = true)
      } text ("Skip synteny alignment and directly merge connected components of BRHs as syntenic anchors.")
      opt[File]("brhs") action { (x, c) =>
        c.copy(brh = x)
      } text ("If BRHs file already exists, provide it here to skip pairwise-ORF alignments.")
      note("\nOPTIONAL: ALIGNMENT PARAMETERS\n")
      opt[Int]('t', "threads") action { (x, c) =>
        c.copy(maxThreads = x)
      } text ("Maximum number of threads to use (default is 1).")
      opt[Int]("kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Kmer sizer used during alignment (default is 11).")
      opt[Int]("minimizer-window") action { (x, c) =>
        c.copy(minimizerWindow = x)
      } text ("Minimizer window size (default is 3).")
      opt[Unit]("dump") action { (x, c) =>
        c.copy(dump = true)
      } text ("Will write intermediate tables to disk (turned off by default).")
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      } text ("Additional log messages.")

    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.database)
      syntenicAnchorsC(config)
    }
  }

  def syntenicAnchorsC(config: Config) = {
    //get all files in database
    val database = config.database.listFiles()
    //fetch all hashtables and repeat expansions:
    //  genome+contig -> Sequence of ORF ids
    val path_hashmap_Z = database.find(_.getName == "global_z.txt").get
    //  ORF ids -> genome+contig
    val path_hashmap_Z_prime = database.find(_.getName == "global_z_prime.txt").get
    // genome+contig -> genome
    val path_hashmap_Y_prime = database.find(_.getName == "global_y_prime.txt").get
    //get potential repeat expansions
    val path_repeats = database.find(_.getName == "repetative_regions.txt").get
    //create output file sof syntenic scores
    //val pw_syntenic_scores = new PrintWriter(config.outputDir + "/syntenic_scores_distribution.txt")
    //char index
    val char_alphabet = 'A' to 'z'

    //Construct y prime hashmap
    val hashmap_Y_prime = {
      openFileWithIterator(path_hashmap_Y_prime).foldLeft(HashMap[String, String]())((map, entry) => {
        val split = entry.split("\t")
        map + (split.head -> split(1))
      })
    }

    //get fasta files of orf sequences
    val sequences = config.database.listFiles().filter(_.getName.endsWith("orfs.sequences.fasta")).par
    //set maximum number of threads to use
    sequences.tasksupport = new ForkJoinTaskSupport(new scala.concurrent.forkjoin.ForkJoinPool(config.maxThreads))
    println(timeStamp + "Found " + sequences.size + " genomes in database")

    //Obtain hashmap H, brhs for a given ORF
    //  TO DO: Post-process BRHs to assure unique sequence mapping per ORF (e.g. an ORF does not have a BRF in
    //  multiple sequences
    val hashmap_H = {
      //for when the user specified BRH file
      if (config.brh != null) {
        println(timeStamp + "User specified BRHs file. Constructing hashmap H from file.")
        verifyFile(config.brh)
        //iterate through file and create hashmap
        openFileWithIterator(config.brh).foldLeft(HashMap[Int, Set[Int]]())((h, line) => {
          val split = line.split("\t")
          h + (split.head.toInt -> split(1).split(",").map(_.toInt).toSet)
        })
      }
      //constructh hashmap H from pairwise ORF alignments
      else {
        println(timeStamp + "Performing " + (sequences.size * sequences.size) + " pairwise ORF-set alignments")
        //obtain H_prime hashmap, best alignments by iterating through them in parallel and merging at the end
        val hashmap_H_prime = {
          //iterate through sequences as query, then as ref for pairwise alignments
          sequences.map(query => {
            progress(1000)
            //collect best alignments for each query sequentially
            sequences.foldLeft(HashMap[Int, Set[Int]]())((h, ref) => {
              if (query == ref) h
              //collect best alignments from pairwise ORF ailgnments
              else {
                collectBestAlignments(ref, query, config.gamma, config.kmerSize, config.minimizerWindow, h)
              }
            }).filterNot(_._2.isEmpty)
            //returns an array of hashmap, combine them in the end
          }).seq.flatten.foldLeft(HashMap[Int, Set[Int]]()) { case (h, (key, value)) => {
            val current = h.getOrElse(key, Set[Int]())
            h + (key -> (current ++ value))
          }
          }
        }

        //dump best alignments
        if (config.dump) {
          val pw = new PrintWriter(config.outputDir + "/best_alignments.txt")
          hashmap_H_prime.foreach(x => pw.println(x._1 + "\t" + x._2.mkString(",")))
          pw.close
        }
        println(timeStamp + "Collected " + hashmap_H_prime.map(_._2.size).sum + " best pairwise ORF alignments")
        println(timeStamp + "Identifying bi-directional alignments (BRHs)")
        //iterate through each query orf and retain only BRHs
        hashmap_H_prime.map { case (query_orf, best_alignments) => {
          //remove cases where there is no bi-direction best alignment (e.g. query <-> ref)
          (query_orf, best_alignments.filter(ref => {
            //attempt to look up current ref's best aligments
            val refs_best_alignments = hashmap_H_prime.get(ref)
            //check whether query is an element of ref's best alignments
            refs_best_alignments != None && refs_best_alignments.get.exists(_ == query_orf)
          }))
        }
        }.filterNot(x => x._2.isEmpty)
      }
    }
    //create brhs file if not provided
    if (config.brh == null) {
      //output file for BRHs; note currently overrides original BRH file
      val pw_brhs = new PrintWriter(config.outputDir + "/brhs.txt")
      hashmap_H.foreach { case (orf, brhs) => pw_brhs.println(Seq(orf, brhs.mkString(",")).mkString("\t")) }
      pw_brhs.close
      //output messages
      println(timeStamp + "Collected " + hashmap_H.map(_._2.size).sum + " BRHs")
    }

    println(timeStamp + "Identifying connected components from BRHs")
    //get connected components of the brhs
    val cc_brhs = getConnectedComponents(hashmap_H, Set(), List()).par
    println(timeStamp + "Found " + cc_brhs.size + " connected components")
    //set max threads to use
    cc_brhs.tasksupport = new ForkJoinTaskSupport(new scala.concurrent.forkjoin.ForkJoinPool(config.maxThreads))

    if (config.dump) {
      val pw_cc = new PrintWriter(config.outputDir + "/ccs.sif")
      cc_brhs.foreach(nodes =>
        nodes.foreach(node => pw_cc.print(node + "\tx\t" + nodes.filterNot(_ == node).mkString("\t") + "\n")))
      pw_cc.close
    }

    println(timeStamp + "Loading hash tables Z, Z', boundaries, and positionals")

    //load hashmaps Z and Z'
    val (boundaries_map, hashmap_Z, hashmap_Z_positional) = {
      //iterate through each entry in the file
      val (b, z, zp) = openFileWithIterator(path_hashmap_Z)
        //and createa a hashmap of lower/upper boundaries based on raw orf IDs
        .foldLeft((HashMap[String, (Int, Int)](),
        //a hashmap with an indexed sequence containing orfs that are elements of that sequence
        HashMap[String, IndexedSeq[Int]](),
        //a hashmap with the orf id as key and it's relative positional index as value
        HashMap[Int, Int]())) { case ((bounds_map, set_map, indexes_map), line) => {
        //split line
        val split = line.split("\t")
        //get orfs ids
        val orfs = split(1).split(",").map(_.toInt)
        //set boundary of the orfs in sequence
        (bounds_map + (split.head -> (orfs.head, orfs.last)),
          //set indexed sequence of all orfs in sequence
          set_map + (split.head -> orfs.toIndexedSeq),
          //set positional index of each ORF in sequence
          orfs.foldLeft((0, indexes_map)) { case ((index, map), orf) => (index + 1, map + (orf -> index)) }._2
        )
      }
      }
      (b.mapValues(x => (x._1 - x._1, x._2 - x._1)), z, zp)
    }

    //load hashmap Z'
    val hashmap_Z_prime =
    //iterate through each entry in file
      openFileWithIterator(path_hashmap_Z_prime).foldLeft(HashMap[Int, String]())((map, line) => {
        //split line
        val split = line.split("\t")
        //set orf key in sequence it belongs to
        map + (split.head.toInt -> split(1))
      })

    //get repetative regions as a map ORF ID -> (Min,Max) of repeat region
    val rer = openFileWithIterator(path_repeats).foldLeft(HashMap[Int, Seq[Int]]())((map, line) => {
      val repeat_region = line.split(",").map(_.toInt).toSeq
      repeat_region.foldLeft(map)((local_map, orf) => local_map + (orf -> repeat_region))
    })

    //log messages
    if (config.verbose) {
      println
      println(timeStamp + "Boundaries map:")
      boundaries_map.foreach(x => println(timeStamp + "--" + x))
      println(timeStamp + "Hashmap Z:")
      hashmap_Z.foreach(x => println(timeStamp + "--" + x))
    }

    /**
      * Function that creates a String of chars for aligning syntenic vectors given a sequence of ORF ids and
      * a map scheme where the key is the orf id and the value is a char
      *
      * @return String
      */
    def constructCharString: (Seq[Int], Map[Int, Char]) => String = (sv, char_map) => sv.map(char_map(_)).mkString

    /**
      * Creates a mapping scheme for a given sequence of ORF ids. For each orf id, obtain the BRHs along with itself
      * and map all BRHs to a common index.
      * NOTE: ASSUMES ALL BRHS SETS AGREE
      *
      * @return Map[Int, Char]
      */
    def constructCharMapping: Seq[Int] => Map[Int, Char] = orfs => {
      orfs.map(orf_id => hashmap_H.getOrElse(orf_id, Set[Int]()) + (orf_id)).toSet.foldLeft((Map[Int, Char](), 0)) {
        case ((cmap, char_index), brh_set) => {
          (brh_set.foldLeft(cmap)((local_cmap, orf_id) => local_cmap + (orf_id -> char_alphabet(char_index))), char_index + 1)
        }
      }._1
    }

    /**
      * Compute synteny scores for a given set of syntenic vectors. Input is a reference generic sytenic vector along
      * with it's inverse complement as a tuple and the target generic sytenic vector in same format
      *
      * @return Int
      */
    def computeSVscore: (Seq[Int], (Seq[Int], Seq[Int])) => Int = (subj, target) => {
      //construct a unified alphabet character map
      val unified_map = constructCharMapping(subj ++ target._1)
      //same as above but for reverse
      val unified_map_reverse = constructCharMapping(subj ++ target._2)
      //create the alphabet string
      val char_string = (constructCharString(subj, unified_map), constructCharString(target._1, unified_map))
      //reverse alphabeet string
      val char_string_reverse = (constructCharString(subj, unified_map_reverse), constructCharString(target._2, unified_map_reverse))
      //println(char_string)
      //println(char_string_reverse)
      //compute edit distances
      val score = levenshteinDistance(char_string._1, char_string._2)
      val score_reverse = levenshteinDistance(char_string_reverse._1, char_string_reverse._2)
      //output smallest of two
      min(score, score_reverse)
    }

    /**
      * Function to obtain the syntenic vectors for a given ORF ID.
      *
      * @return
      */
    def getSyntenicVector: Int => ((Seq[Int], Seq[Int]), (Seq[Int], Seq[Int]), (Seq[Int], Seq[Int])) = orf_id => {
      //obtain native sequence/chrm of given ORF
      val native_sequence = hashmap_Z_prime(orf_id)
      //obtain the entire sequence of ORFs in the native sequence
      val sequence_orfs = hashmap_Z(native_sequence)
      //index of (potentially) computational start/end ORFs
      val boundaries = boundaries_map(native_sequence)
      //check to see if orf is part of a repeat expansion
      val repeat = rer.get(orf_id)
      //get index of the given ORF in the seuence of ORFs. in the case of a repeat expansion, will get left and
      // right-most indeces
      val (left_most_index, right_most_index) = {
        //if orf id is not part of repeat expansion, get the index, right index is -1 (i.e. ignore)
        if(repeat.isEmpty) {
          //get index
          val index = hashmap_Z_positional(orf_id)
          //return both
          (index, index)
        }
        //orf is is part of repeat expansion, get left and right-most indeces
        else (hashmap_Z_positional(repeat.get.head), hashmap_Z_positional(repeat.get.last))
      }

      //function to determine whether current orf is at the ends of a sequence based on flanking window parameter
      def atEnds(): Boolean =
        config.flankingWindow - left_most_index < 0 || config.flankingWindow + right_most_index > boundaries._2

      //log messages
      if (config.verbose) {
        println
        println(timeStamp + "ORF ID: " + orf_id)
        println(timeStamp + "Native sequence: " + native_sequence)
        println(timeStamp + "ORF sequence: " + sequence_orfs)
        println(timeStamp + "Boundaries: " + boundaries)
        println(timeStamp + "ORF index(s): " + (left_most_index, right_most_index))
      }

      /**
        * Function to construct the left and right syntenic vectors
        *
        * @return
        */
      def constructFlankingSV(): (Seq[Int], Seq[Int], Seq[Int]) = {
        //iterate from 1 to flanking window and add accordingly
        val (lsv, rsv) = (1 to config.flankingWindow).foldLeft((Seq[Int](), Seq[Int]()))((sv, offset) => {
          //update left sv
          val left = {
            //compute different between current index and current offset
            val left_index = (left_most_index - offset)
            //far enough from the edge, add accordingly
            if (left_index >= 0) Seq(sequence_orfs(left_index)) ++ sv._1
            //need to circularize and add accordingly
            else Seq(sequence_orfs(boundaries._2 + (left_index + 1))) ++ sv._1
          }
          //update right sv
          val right = {
            //same as above
            val right_index = (right_most_index + offset)
            //same as above
            if (boundaries._2 - right_index >= 0) sv._2.:+(sequence_orfs(right_index))
            //same as above
            else sv._2.:+(sequence_orfs(boundaries._1 + ((boundaries._2 - right_index + 1) * -1)))
          }
          (left, right)
        })
        //return left, general, and right
        (lsv, lsv.takeRight(config.flankingWindow / 2).:+(orf_id) ++ rsv.take(config.flankingWindow / 2), rsv)
      }

      /**
        * Function to construct general syntenic vector for non-circular genomes
        *
        * @return Seq[Int]
        */
      def constructGeneralSV(): Seq[Int] = {
        //calculate the left most offset and potential carryon offset for the right
        val (left_most_offset, right_carryon) = {
          //the difference between the flanking window and the current index
          val offset = left_most_index - config.flankingWindow / 2
          //too close to the edge, add left-most reachable index the rest as carry on for the right
          if (offset < 0) (left_most_index - boundaries._1, offset * -1)
          //far enough from the edge, add accordinlgy
          else (config.flankingWindow / 2, 0)
        }
        //same as above but for the right
        val (right_most_offset, left_carryon) = {
          val offset = boundaries._2 - (right_most_index + config.flankingWindow / 2)
          if (offset < 0) (boundaries._2 - right_most_index, offset * -1)
          else (config.flankingWindow / 2, 0)
        }
        //functional way to construct general sv:
        //  iterate through left-most offset to 1 and add to sequence
        (1 to (left_most_offset + left_carryon)).reverse.foldLeft(Seq[Int]())((left_sv, offset) => {
          left_sv.:+(sequence_orfs(left_most_index - offset))
        }) ++
          // iterate from 1 to right-most offset and add to sequence
          (1 to (right_most_offset + right_carryon)).foldLeft(Seq[Int]())((right_sv, offset) => {
            right_sv.:+(sequence_orfs(right_most_index + offset))
          })
      }

      //when the genome is circular
      if (config.isCircular) {
        //compute syntenic vectors
        val (left_sv, general_sv, right_sv) = constructFlankingSV()
        //return left general and right along with their corresponding reverse
        ((left_sv, right_sv.reverse), (general_sv, general_sv.reverse), (right_sv, left_sv.reverse))
      }
      else {
        //check if the current ORF is too close to the ends, therefore force general window only
        if (atEnds()) {
          val general_sv = constructGeneralSV()
          val (left_sv, right_sv) = (Seq[Int](), Seq[Int]())
          ((left_sv, right_sv.reverse), (general_sv, general_sv.reverse), (right_sv, left_sv.reverse))
        }
        //create svs as normal
        else {
          val (left_sv, general_sv, right_sv) = constructFlankingSV()
          ((left_sv, right_sv.reverse), (general_sv, general_sv.reverse), (right_sv, left_sv.reverse))
        }
      }
    }

    //create output file
    val pw_syntenic_anchors = new PrintWriter(config.outputDir + "/syntenic_anchors.txt")

    println(timeStamp + "Computing syntenic anchors in each connected component")
    //iterate through each connected component
    cc_brhs.foreach(cc => {
      //automatically merge ccs into a syntenic anchor if parameter is turned on
      if (config.mergeCC) cc.foreach(x => pw_syntenic_anchors.println(x + "\t" + cc.filterNot(_ == x).mkString(",")))
      //else proceed to synteny alignment
      else {
        //group by genome
        val cc_grouped = cc.groupBy(x => hashmap_Y_prime(hashmap_Z_prime(x))).toList
          //iterate through each genome and set of ORFs
          //.mapValues(group => {
          //check if group is part of repeat expansion. ASSUMES THAT IF ONE ORF IS AN RER, THAN SO ARE REST
          //val isRepeatExpansion = group.exists(rer.contains(_))
          //if it's a repeat expansion, arbitraily get on
          //if (isRepeatExpansion) group.take(1) else group
        //}).toList
        //all brhs are unique (no multiple brhs in the same genome)
        if(cc_grouped.forall(_._2.size == 1))
          cc.foreach(x => pw_syntenic_anchors.println(x + "\t" + cc.filterNot(_ == x).mkString(",")))
        else {
        println(timeStamp + "--Found multiple candidate RBHs. Performing synteny alignment.")
          println(cc)
        //pairwise comparison of all orfs
        val syntenic_anchors = cc.foldLeft(List[((Int, Int), Int)]())((synteny_scores, subject_orf) => {
          //construct syntenic vectors for current orf
          val (subject_left_sv, subject_general_sv, subject_right_sv) =  getSyntenicVector(subject_orf)
          //iterate through target orfs
          cc.foldLeft(synteny_scores)((local_synteny_scores, target_orf) => {
            //skip if orfs are from same genome
            if (hashmap_Y_prime(hashmap_Z_prime(subject_orf)) == hashmap_Y_prime(hashmap_Z_prime((target_orf)))) local_synteny_scores
            //calculate synteny score
            else {

              //println(subject_orf, target_orf)
              //println((subject_left_sv, subject_general_sv, subject_right_sv))
              //get syntenic vectors of target
              val (target_left_sv, target_general_sv, target_right_sv) = getSyntenicVector(target_orf)
              //println((target_left_sv, target_general_sv, target_right_sv))
              //get smallest synteny score
              val synteny_score = List(computeSVscore(subject_left_sv._1, target_left_sv),
                computeSVscore(subject_general_sv._1, target_general_sv),
                computeSVscore(subject_right_sv._1, target_right_sv)).min
              //update synteny score
              local_synteny_scores.:+((subject_orf, target_orf), synteny_score)
            }
          })
        })
          //process connected components into set of syntenic anchors
          //remove high-distant orfs
          .filter(_._2 <= 1)
          //group by subject orf
          .groupBy(_._1._1)
          //group values by genome and in the case of multiple orfs for same genome, retain lowest scoring one
          .mapValues(_.groupBy(x => hashmap_Y_prime(hashmap_Z_prime(x._1._2))).mapValues(_.minBy(_._2)._1._2))
          //merge all orfs into one set, creating the syntenic anchors
          //--to do: see if this needs partitioning identify highly connected subgraphs?
          .foldLeft(Set[Int]())((b, a) => a._2.foldLeft(b)((c, d) => c + d._2) + a._1)
        //iterate and output to file
        syntenic_anchors.foreach(x =>
          pw_syntenic_anchors.println(Seq(x, syntenic_anchors.filterNot(_ == x).mkString(",")).mkString("\t")))
      }}
    })
    pw_syntenic_anchors.close
    println(timeStamp + "Succesfully completed!")
  }

  /**
    * Function to obtain the name of a file given a File object
    *
    * @return String
    */
  def getName: File => String = name => name.getName.substring(0, name.getName.lastIndexOf("."))

  /**
    * Function to compute the mean displacement
    *
    * @return
    */
  def mean: Seq[Double] => Double = displacements => displacements.sum / displacements.size

  /**
    * Function to compute the levenshtein distance. Shout-out to:
    *
    * @param a
    * @param b
    * @tparam A
    * @return
    */
  def levenshteinDistance[A](a: Iterable[A], b: Iterable[A]) = {
    ((0 to b.size).toList /: a) ((prev, x) =>
      (prev zip prev.tail zip b).scanLeft(prev.head + 1) {
        case (h, ((d, v), y)) => min(min(h + 1, v + 1), d + (if (x == y) 0 else 1))
      }) last
  }

}
