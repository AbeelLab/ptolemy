package syntenic_anchors

import java.io.{File, PrintWriter}

import utilities.FileHandling.timeStamp
import utilities.FileHandling.{openFileWithIterator, tLines, verifyDirectory, verifyFile}
import utilities.GFFutils
import utilities.MinimapUtils

import scala.collection.immutable.{HashMap, HashSet}
import atk.Tool.progress

import scala.annotation.tailrec


/**
  * Author: Alex N. Salazar
  * Created on 14-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SyntenicAnchors extends tLines with GFFutils with MinimapUtils {

  case class Config(
                     database: File = null,
                     outputDir: File = null,
                     gamma: Double = 0.75,
                     f: Int = 10,
                     verbose: Boolean = false,
                     dump: Boolean = false,
                     syntenicFraction: Double = 0.5,
                     kmerSize: Int = 11,
                     minimizerWindow: Int = 5,
                     brh: File = null,
                     debug: Boolean = false
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("syntenic-anchors") {
      opt[File]("db") required() action { (x, c) =>
        c.copy(database = x)
      } text ("Directory path of database (e.g. output directory of 'extract' module).")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory.")
      note("\nOPTIONAL\n")
      opt[File]("brhs") action { (x, c) =>
        c.copy(brh = x)
      } text ("If best reciprocal hits file already exists, provide it here to skip pairwise-ORF alignments.")
      opt[Double]("gamma") action { (x, c) =>
        c.copy(gamma = x)
      } text ("Gamma-value for minimum coverage of ORF alignment (default is 0.3).")
      opt[Int]("f") action { (x, c) =>
        c.copy(f = x)
      } text ("Flanking window (default is 10).")
      opt[Double]("beta") action { (x, c) =>
        c.copy(syntenicFraction = x)
      } text ("Maximum allowed (default is 0.5).")
      opt[Int]("kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Kmer sizer used during alignment (default is 11).")
      opt[Int]("minimizer-window") action { (x, c) =>
        c.copy(minimizerWindow = x)
      } text ("Minimizer window size (default is 5).")
      opt[Unit]("dump") action { (x, c) =>
        c.copy(dump = true)
      } text ("Will write intermediate tables to disk (turned off by default).")
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      }

    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.database)
      getBRHs(config)
    }
  }

  def getBRHs(config: Config) = {
    val available_corse = Runtime.getRuntime().availableProcessors()
    val database = config.database.listFiles()
    val path_hashmap_Z = database.find(_.getName == "global_z.txt").get
    val path_hashmap_Z_prime = database.find(_.getName == "global_z_prime.txt").get
    val path_hashmap_Y_prime = database.find(_.getName == "global_y_prime.txt").get
    val path_repeats = database.find(_.getName == "repetative_regions.txt").get
    val pw_syntenic_scores = new PrintWriter(config.outputDir + "/syntenic_scores_distribution.txt")
    //construct y prime hashmap
    val hashmap_Y_prime = {
      openFileWithIterator(path_hashmap_Y_prime).foldLeft(HashMap[String, String]())((map, entry) => {
        val split = entry.split("\t")
        map + (split.head -> split(1))
      })
    }
    //get fasta files of orf sequences
    val sequences = config.database.listFiles().filter(_.getName.endsWith("fasta"))
    println(timeStamp + "Found " + sequences.size + " genomes in database")
    //obtain hashmap H, brhs for a given ORF
    val hashmap_H = {
      //for when the user specified BRh file
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
        println(timeStamp + "Performing " + (sequences.size * sequences.size) + " pairwise ORF-set alignments (in parallel" +
          " if multiple cores are available): ")
        //obtain H_prime hashmap, best alignments by iterating through them in parallel and merging at the end
        val hashmap_H_prime = {
          //iterate through sequences as query, then as ref for pairwise alignments
          sequences.par.map(query => {
            //collect best alignments for each query sequentially
            sequences.foldLeft(HashMap[Int, Set[Int]]())((h, ref) => {
              if(available_corse < 10 ) progress(100) else progress (1000)
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
          hashmap_H_prime.foreach(x => pw.println(Seq(x._1, x._2.mkString(",")).mkString("\t")))
          pw.close
        }
        println(timeStamp + "Collected " + hashmap_H_prime.map(_._2.size).sum + " best pairwise ORF alignments")
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
    val pw_brhs = new PrintWriter(config.outputDir + "/brhs.txt")
    hashmap_H.foreach { case (orf, brhs) => pw_brhs.println(Seq(orf, brhs.mkString(",")).mkString("\t")) }
    pw_brhs.close
    println(timeStamp + "Collected " + hashmap_H.map(_._2.size).sum + " BRHs")
    println(timeStamp + "Loading hash tables Z, Z', boundaries, and positionals")
    //load hashmaps Z and Z'
    val (boundaries_map, hashmap_Z, hashmap_Z_positional) =
    //iterate through each entry in the file
      openFileWithIterator(path_hashmap_Z)
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

    /**
      * Function to check whether an orf is part of a repetative region. If so, will return left/right most ORF id of
      * the region. If not, returns a tuple of itself
      *
      * @return (Int,Int)
      */
    def getRepeatRegionsEnds: Int => (Int, Int) = orf_id => {
      val get = rer.get(orf_id)
      if (get == None) (orf_id, orf_id) else (get.get.head, get.get.last)
    }

    /**
      * Function that tries to give the corresponding orfs of the same rank if either of the two given orf IDs exist
      * within a repetative region. Does this by trying to get the query orf that is in the same rank (index) as that
      * of the reference orf. Returns a tuple with the first being a boolean if both orfs can be properly aligned
      * with their ranks.
      *
      * @return
      */
    def getRepeatRank: (Int, Int) => (Boolean, Int, Int) = (r_orf, q_orf) => {
      //get r_orf's repetative region, if it exists
      val r_orf_repeat_region = rer.getOrElse(r_orf, Seq(r_orf))
      //get q_orf's repetative region, if it exists
      val q_orf_repeat_region = rer.getOrElse(q_orf, Seq(q_orf))
      //get the rank of the reference orf
      val r_orf_rank = r_orf_repeat_region.indexOf(r_orf)
      //if their ranks can't be properly match, return accordingly
      if (r_orf_rank + 1 > q_orf_repeat_region.size) (false, r_orf, q_orf)
      else (true, r_orf, q_orf_repeat_region(r_orf_rank))
    }

    /**
      * Function to compute positional displacement given two ORF IDs.
      * Note: using raw ORF ids as proxies of positional indexes
      *
      * @return Positional displacement of two ORFs
      */
    def computeDisplacement: (Int, Int) => Double = (orf_a, orf_b) => {
      val d = Math.abs(hashmap_Z_positional(orf_a) - hashmap_Z_positional(orf_b))
      if (config.verbose) println("------Displacement for " + (orf_a, orf_b) + " is " + d)
      d
    }

    /**
      * Function to compute the positional displacement of the BRH of orf_a in respects to some ORF, orf_b.
      * Controls cases when multiple or no local BRHs exist within a flanking window of size f.
      *
      * @return Positional displacement of the BRH of orf_a in respects to orf_b
      */
    def phi: (Int, Int, Boolean) => Double = (orf_a, orf_b, within_window) => {
      //get corresponding sequence of orf_b
      val orf_b_sequence = hashmap_Z_prime(orf_b)
      //obtain all BRHs within a flanking window f
      val L = {
        //try to get brhs for current orf
        val brhs = hashmap_H.get(orf_a)
        //if none exists (e.g. unique orf) return empty
        if (brhs == None) List[Int]()
        else brhs.get.filter(x => {
          //same sequence as orf_b
          hashmap_Z_prime(x) == orf_b_sequence &&
            //and within f position from orf_b if desired; note using raw orf ids as positional proxies
            {
              if (within_window) Math.abs(orf_b - x) <= config.f else true
            }
        })
      }
      //when L is empty, return BRH penalty, else compute minimum of positional displacement
      val d = L.size match {
        case 0 => -(config.f*config.syntenicFraction)
        case _ => L.map(orf_l => computeDisplacement(orf_l, orf_b)).min
      }
      if (config.verbose) println("------Phi value is " + d + " with L size of " + L.size + " for orf " + orf_a)
      d
    }

    /**
      * Function to get positional index boundaries for a given orf id
      *
      * @return (Lower,Upper) boundaries
      */
    def getBoundaries: Int => (Int, Int) = orf_id => boundaries_map(hashmap_Z_prime(orf_id))

    /**
      * Generic function to determine if a given value is contained withing a given range
      *
      * @return Boolean
      */
    def isContained: (Int, (Int, Int)) => Boolean = (value, range) => range._1 <= value && range._2 >= value

    /**
      * Function to determine whether an orf will be within boundaries x positions away
      * Note: x can be negative
      *
      * @return
      */
    def isWithinBoundaries: (Int, Int) => Boolean = (orf_id, x) => {
      val statement = isContained(orf_id + x, getBoundaries(orf_id))
      if (config.verbose && !statement) println("------ORF " + orf_id + " not within bounds for " + x)
      statement
    }

    /**
      * Function to determine the left and right most positions for a ORF that is near the ends of a sequence
      *
      * @return
      */
    def getLeftRightMostReach: Int => (Int, Int) = orf_id => {
      //get general window size
      val general_window_size = if (config.f % 2 == 0) config.f / 2 else (config.f / 2) + 1
      //get contextual position of orf
      val position = hashmap_Z_positional(orf_id)
      //get boundaries for the sequence of orf
      val boundaries = getBoundaries(orf_id)
      //get left and right offsets
      val left_offset = position - boundaries._1
      val right_offset = boundaries._2 - position
      //if both are reachable, return
      if (left_offset >= 0 && right_offset >= 0) (general_window_size, general_window_size)
      else {
        //assuming one of the two sides is reachable, throw warning if this is not true
        if (right_offset < 0 && left_offset < 0) println("WARNING: could not find general offset for ORF " + orf_id)
        //add offset accordingly to one of the sides
        val final_left_offset = if (right_offset < 0) general_window_size + (-right_offset) else general_window_size
        val final_right_offset = if (left_offset < 0) general_window_size + (-left_offset) else general_window_size
        (final_left_offset, final_right_offset)
      }
    }

    def getEndMostReach: (Int, Int, Int) => (Int,Int) = (orf_a, orf_b, direction) => {
      //get general window size
      val general_window_size = if (config.f % 2 == 0) config.f / 2 else (config.f / 2) + 1
      //get contextual position of orf
      val (position_a, position_b) = (hashmap_Z_positional(orf_a), hashmap_Z_positional(orf_b))
      //get boundaries for the sequence of orf
      val (boundaries_a, boundaries_b) = (getBoundaries(orf_a), getBoundaries(orf_b))
      val min_offset = {
        direction match {
          case 0 => List((position_a - general_window_size) - hashmap_Z_positional(boundaries_a._1),
            (position_b - general_window_size) - hashmap_Z_positional(boundaries_b._1)).min
          case 1 => List(hashmap_Z_positional(boundaries_a._2) - (position_a + general_window_size),
            hashmap_Z_positional(boundaries_b._2) - (position_b + general_window_size)).min
        }
      }
      if(min_offset >= 0) (general_window_size, 0) else (general_window_size + min_offset, -min_offset)
    }


    /**
      * Function to compute minimum observed positional displacement of a BRH that is x positions away.
      * Note: x can be negative
      *
      * @return Minimum positional displacement of a BRH x positions away
      */
    def psi: (Int, Int, Int, Int, Boolean) => Double = (orf_i, orf_j, x, y, w) => {
      if (config.verbose) println("----Computing displacements with " + x + " value")
      //compute displacement in the perspective or orf i
      val perspective_orf_i = {
        if (!isWithinBoundaries(orf_i, x)) (config.f*config.syntenicFraction)
        else {
          val phi_score = phi(orf_i + x, orf_j, w)
          val displacement_score = computeDisplacement(orf_i + x, orf_i)
          Math.abs(displacement_score - phi_score)
        }
      }
      //compute displacement in the perspective or orf i
      val perspective_orf_j = {
        if (!isWithinBoundaries(orf_j, y)) (config.f*config.syntenicFraction)
        else {
          val phi_score = phi(orf_j + y, orf_i, w)
          val displacment_score = computeDisplacement(orf_j + y, orf_j)
          Math.abs(displacment_score - phi_score)
        }
      }
      if (config.verbose) println("----Found the following values " + (perspective_orf_i, perspective_orf_j) + " " +
        "translating to " + List(perspective_orf_i, perspective_orf_j).max)
      List(perspective_orf_i, perspective_orf_j).max
    }

    def constructSyntenicVectors: (Int, Int) =>
      (Seq[Double], Seq[Double], Seq[Double], Seq[Double], Option[(Int,Int)]) = (orf_a, orf_b) => {
      if (config.verbose) println("--Using following positions away for left SV: " + (0 to config.f - 1))
      //get respective orf ids if they belong to repetative regions
      val (orf_a_left, orf_a_right) = getRepeatRegionsEnds(orf_a)
      val (orf_b_left, orf_b_right) = getRepeatRegionsEnds(orf_b)
      //check if all nodes are reachable with default flanking sizes
      val is_all_reachable =
        List(orf_a_left, orf_a_right, orf_b_left, orf_b_right).forall(x => isWithinBoundaries(x, config.f))
      //get flanking window sized depending on distance to the ends
      val (left_flank_size, right_flank_size, general_flank) = {
        //all nodes are reachable with default window sizes
        if (is_all_reachable) (config.f, config.f, None)
          //need to use general flanking window mode
        else {
          if (config.verbose) println("--Turning on general flank mode")
          //get the minimum left-most position and the offset, if any
          val (left_window, left_offset) = getEndMostReach(orf_a_left, orf_b_left, 0)
          //get the minimum right-most position and the offset, if any
          val (right_window, right_offset) = getEndMostReach(orf_a_right, orf_b_right, 1)
          if(config.verbose) println("--Found the following parameters: " + (left_window, left_offset, right_window,
            right_offset))
          //construct general window
          val (left, right) = (left_window + right_offset, right_window + left_offset)
          (left, right, Some((left, right)))
        }
      }
      //compute left syntenic vector
      val (sv_l, sv_l_prime) = (0 to (left_flank_size - 1)).foldLeft((Seq[Double](), Seq[Double]())) { case ((sv, svi), p) => {
        val (x, y) = (-(config.f - p), -(config.f - p))
        (sv.:+(psi(orf_a_left, orf_b_left, x, y, true)), svi.:+(psi(orf_a_left, orf_b_left, x, -y, false)))
      }
      }
      if (config.verbose) println("--Using following positions away for right SV: " + (1 to config.f))
      //compute left syntenic vector
      val (sv_r, sv_r_prime) = (1 to right_flank_size).foldLeft((Seq[Double](), Seq[Double]())) { case ((sv, svi), p) => {
        val (x, y) = (p, p)
        (sv.:+(psi(orf_a_right, orf_b_right, x, y, true)), svi.:+(psi(orf_a_right, orf_b_right, x, -y, false)))
      }
      }
      if (config.verbose) {
        if (general_flank == None){
          println("--Syntenic vectors: " + (sv_l, sv_r))
        println("--Inverse syntenic vectors: " + (sv_l_prime, sv_r_prime))
      } else {
        println("--General flank: " + general_flank.get + "\t" + (sv_l ++ sv_r))
        println("--Inverse general flank " + general_flank.get + "\t" + (sv_l_prime, sv_r_prime))
      }
    }
      (sv_l, sv_r, sv_l_prime, sv_r_prime, general_flank)
    }

    println(timeStamp + "Computing syntenic scores for " + hashmap_H.size + " BRH sets")
    val syntenic_anchor_pairs = hashmap_H.toList.par.flatMap { case (orf, brhs) => {
      progress(10000)
      /**
        * Method that iterate through a list of BRHs and retains the best syntenic score along with corresponding ORF
        *
        * @param _brhs
        * @param min_score
        * @param min_orfs
        * @return
        */
      @tailrec def computeSyntenicScore(_brhs: Set[Int],
                               min_score: Option[Double],
                               min_orfs: List[(Int, Seq[Double])]): (Option[Double], List[(Int, Seq[Double])]) = {
        if (_brhs.isEmpty) (min_score, min_orfs)
        else {
          if (config.verbose) println("Scoring BRH: " + (orf, _brhs.head))
          //get the corresponding ranked ORF if either belong to a repetative region
          val (is_same_rank, ranked_orf, ranked_brh) = getRepeatRank(orf, _brhs.head)
          //if they are not the same rank, move on
          if (!is_same_rank) computeSyntenicScore(_brhs.tail, min_score, min_orfs)
          else {
            //compute left and right syntenic vectors
            val (sv_l, sv_r, sv_l_prime, sv_r_prime, general_flank) = constructSyntenicVectors(ranked_orf, ranked_brh)
              //compute syntenic scores
              val (min_syntenic_score, all_scores) = {
                if (general_flank == None) {
                  //compute synteny scores
                  val (sv_l_score, sv_r_score, sv_l_prime_score, sv_r_prime_score) =
                    (mean(sv_l), mean(sv_r), mean(sv_l_prime), mean(sv_r_prime))
                  //output to file
                  if(config.dump) pw_syntenic_scores.println(Seq(ranked_orf, ranked_brh.toString, sv_l_score,
                    sv_r_score, sv_l_prime_score, sv_r_prime_score).mkString("\t"))
                  //create scores for all vectors
                  val all_scores = Seq(sv_l_score, sv_r_score, sv_l_prime_score, sv_r_prime_score)
                  //return
                  (Math.pow(all_scores.min, 1.5), all_scores)
                }
                else {
                  //compute syntenic scores
                  val (sv_general_score) = mean(sv_l ++ sv_r)
                  if(config.verbose) println("--General score: " + sv_general_score)
                  if (config.dump) pw_syntenic_scores.println(
                    Seq(ranked_orf, ranked_brh.toString, sv_general_score, "f:" + general_flank.get).mkString("\t"))
                  (Math.pow(Seq(sv_general_score).min, 1.5), Seq(sv_general_score))
                }
              }
              if (config.verbose) println("--Syntenic score is: " + min_syntenic_score)
              //update accordingly if the syntenic score is lower than previously observed
              if (min_score == None || min_syntenic_score < min_score.get) {
                computeSyntenicScore(_brhs.tail, Option(min_syntenic_score), List((ranked_brh, all_scores)))
              } else if(min_syntenic_score == min_score.get){
                computeSyntenicScore(_brhs.tail, Option(min_syntenic_score), (ranked_brh, all_scores) :: min_orfs)
              }
              else computeSyntenicScore(_brhs.tail, min_score, min_orfs)
          }
        }
      }
      //retain syntenic anchor per genome
        brhs.groupBy(x => hashmap_Y_prime(hashmap_Z_prime(x))).foldLeft(List[(Int,Int)]()){case(sa, local_brhs) => {
        //fetch best syntenic score and the corresponding orf
        val (best_syntenic_score, best_orfs) = computeSyntenicScore(local_brhs._2, None, List())
        //control when there are unique sequences and when the syntenic scores is above specified threshold
        if (best_syntenic_score == None) sa
        else if (!(best_syntenic_score.get <= (config.f * config.syntenicFraction))) sa
        else {
          //check the number of orfs
          best_orfs.size match {
              //only one orf with the best score
            case 1 =>{
              if (config.verbose) println("Classified as a syntenic anchor: " + (local_brhs._1, (orf, best_orfs.head._1)))
              //(local_brhs._1, (orf, best_orfs.head._1)) :: sa
              (orf, best_orfs.head._1) :: sa
            }
              //multiple orfs with the same best score
            case _ => {
              if(config.verbose) println("Multiple ORFs with the best syntenic scores: " + (local_brhs._1, (orf, best_orfs)))
              //function to the get the scores of the inverses
              def getInverse: Seq[Double] => (Double, Double) = vector => {
                if (vector.size == 1) (vector.head, vector.head)
                else (vector.takeRight(2).take(1).head, vector.takeRight(1).head)
              }
              //rank the multiple best orfs based on: sum of left/right/general flank, left-inverse, right-inverse
              val best_orf = best_orfs.sortBy(x => {
                val inverses = getInverse(x._2)
                (x._2.take(2).sum, inverses._1, inverses._2)
              }).head
              //(local_brhs._1, (orf, best_orf._1)) :: sa
              (orf, best_orf._1) :: sa
            }
          }
        }
      }}
    }
    }
    //dump syntenic pairs
    if (config.dump) {
      val pw = new PrintWriter(config.outputDir + "/syntenic_pairs.txt")
      syntenic_anchor_pairs.foreach(x => pw.println(Seq(x._1, x._2).mkString("\t")))
      pw.close
    }
    //get set of orfs that are part of the syntenic anchors and compute size of the set
    val total_orfs_in_syntenic_anchor_pairs = syntenic_anchor_pairs.seq.flatMap(x => List(x._1, x._2)).toSet
    println(timeStamp + "Found " + total_orfs_in_syntenic_anchor_pairs.size + " involved in " + syntenic_anchor_pairs.size +
      " syntenic anchor pairs")
    println(timeStamp + "Constructing syntenic anchor graph.")
    //construct syntenic anchor graph
    val syntenic_anchor_graph =
      syntenic_anchor_pairs.foldLeft(HashMap[Int, Set[Int]]()) { case (map_syntenic_anchors, (orf_a, orf_b)) => {
        val current_a = map_syntenic_anchors.getOrElse(orf_a, Set[Int]())
        val current_b = map_syntenic_anchors.getOrElse(orf_b, Set[Int]())
        val update1 = (map_syntenic_anchors + (orf_a -> current_a.+(orf_b)))
        update1 + (orf_b -> current_b.+(orf_a))
      }
      }

    /**
      * Method to find connectect components in the syntenic anchor graph using breadth-first search traversal
      *
      * @param nodes
      * @param ccs
      * @return
      */
    @tailrec def identifySyntenicAnchorCC(nodes: Set[Int], ccs: List[Set[Int]]): List[Set[Int]] = {
      /**
        * Method to perform breadth-first search on the syntenic-anchor graph
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
            val neighbours = syntenic_anchor_graph.getOrElse(head, Set[Int]()).toSeq.filterNot(localVisited.contains(_))
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
        identifySyntenicAnchorCC(nodes.tail -- cc, cc :: ccs)
      }
    }

    //create output file
    val pw_syntenic_anchors = new PrintWriter(config.outputDir + "/syntenic_anchors.txt")
    //get connected components in the syntenic anchor graph
    val syntenic_anchors = identifySyntenicAnchorCC(total_orfs_in_syntenic_anchor_pairs, List())
    println(timeStamp + "Found " + syntenic_anchors.size + " connected components.")
    println(timeStamp + "Writing to disk.")
    syntenic_anchors.foreach(cc => {
      cc.foreach(node => {
        pw_syntenic_anchors.println(Seq(node, cc.filterNot(_ == node).mkString(",")).mkString("\t"))
      })
    })
    pw_syntenic_anchors.close
    pw_syntenic_scores.close
    println(timeStamp + "Successfully completed!")

  }

  def getName: File => String = name => name.getName.substring(0, name.getName.lastIndexOf("."))

  def mean: Seq[Double] => Double = displacements => displacements.sum / displacements.size
}
