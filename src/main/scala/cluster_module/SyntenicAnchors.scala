package cluster_module

import java.io.{File, PrintWriter}

import utilities.FileHandling.timeStamp
import utilities.FileHandling.{tLines, verifyDirectory, openFileWithIterator}
import utilities.GFFutils
import utilities.ORFalignments
import utilities.BRHs

import scala.collection.immutable.HashMap
import atk.Tool.progress


/**
  * Author: Alex N. Salazar
  * Created on 14-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SyntenicAnchors extends tLines with GFFutils with ORFalignments with BRHs {

  case class Config(
                     database: File = null,
                     outputDir: File = null,
                     gamma: Double = 0.75,
                     dpenalty: Double = 0.5,
                     f: Int = 10,
                     alpha: Int = 2,
                     verbose: Boolean = false,
                     dump: Boolean = false,
                     syntenicFraction: Double = 0.75,
                     kmerSize: Int = 11,
                     minimizerWindow: Int = 1,
                     debug: Boolean = false
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("get-brhs") {
      opt[File]("db") required() action { (x, c) =>
        c.copy(database = x)
      } text ("Directory path of database (e.g. output directory of 'extract' module).")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory.")
      note("\nOPTIONAL\n")
      opt[Double]("gamma") action { (x, c) =>
        c.copy(gamma = x)
      } text ("Gamma-value for minimum coverage of ORF alignment (default is 0.75).")
      opt[Int]("flanking-window") action { (x, c) =>
        c.copy(f = x)
      } text ("Flanking window (default is 10).")
      opt[Int]("alpha") action { (x, c) =>
        c.copy(alpha = x)
      } text ("BRH penalty (default to 2).")
      opt[Double]("displacement-penalty") action { (x, c) =>
        c.copy(dpenalty = x)
      } text ("Displacement penalty (default to 0.5).")
      opt[Double]("syntenic-fraction") action { (x, c) =>
        c.copy(syntenicFraction = x)
      } text ("Maximum allowed (default to 0.75).")
      opt[Int]("kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Kmer sizer used during alignment (default is 11).")
      opt[Int]("minimizer-window-size") action { (x, c) =>
        c.copy(minimizerWindow = x)
      } text ("Minimizer window size (default is 1).")
      opt[Unit]("dump") action { (x, c) =>
        c.copy(dump = true)
      } text("Will write intermediate tables to disk (turned off by default).")
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
    val path_hashmap_Z = config.database.listFiles().find(_.getName == "global_z.txt").get
    val path_hashmap_Z_prime = config.database.listFiles().find(_.getName == "global_z_prime.txt").get
    val path_hashmap_Y_prime = config.database.listFiles().find(_.getName == "global_y_prime.txt").get
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
    println("Found " + sequences.size + " genomes in database" + timeStamp)
    println("Performing " + (sequences.size * sequences.size) + " pairwise ORF-set alignments:" + timeStamp)
    //obtain hashmap H, brhs for a given ORF
    val hashmap_H = {
      //obtain H_prime hashmap, best alignments
      val hashmap_H_prime =
      //iterate through sequences as query, then as ref for pairwise alignments
        sequences.foldLeft(HashMap[Int, List[Int]]())((_hashmap_H, query) => {
          sequences.foldLeft(_hashmap_H)((local_hashmap_H, ref) => {
            progress(1)
            if (query == ref) local_hashmap_H
            //collect best alignments from pairwise ORF ailgnments
            else {
              collectBestAlignments(ref, query, config.gamma, config.kmerSize, config.minimizerWindow, local_hashmap_H)
            }
          })
        })
      //dump best alignments
      if(config.dump) {
        val pw = new PrintWriter(config.outputDir + "/best_alignments.txt")
        hashmap_H_prime.foreach(x => pw.println(Seq(x._1, x._2.mkString(",")).mkString("\t")))
        pw.close
      }
      println("Collected " + hashmap_H_prime.map(_._2.size).sum + " best pairwise ORF alignments" + timeStamp)
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
      }
    }
    //dump brhs
    if (config.dump) {
      val pw = new PrintWriter(config.outputDir + "/brhs.txt")
      hashmap_H.foreach { case (orf, brhs) => pw.println(Seq(orf, brhs.mkString(",")).mkString("\t")) }
      pw.close
    }
    println("Collected " + hashmap_H.map(_._2.size).sum + " BRHs" + timeStamp)
    println("Loading hash tables Z, Z', boundaries, and positionals" + timeStamp)
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
            {if(within_window) Math.abs(orf_b - x) <= config.f else true}
        })
      }
      //when L is empty, return BRH penalty, else compute minimum of positional displacement
      val d = L.size match {
        case 0 => config.alpha
        case _ => L.map(orf_l => computeDisplacement(orf_l, orf_b)).min
      }
      if (config.verbose) println("------Phi value is " + d + " with L size of " + L.size + " for orf " + orf_a)
      d
    }

    /**
      * Function to determine whether an orf will be within boundaries x positions away
      * Note: x can be negative
      *
      * @return
      */
    def isWithinBoundaries: (Int, Int) => Boolean = (orf_id, x) => {
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

      val statement = isContained(orf_id + x, getBoundaries(orf_id))
      if (config.verbose && !statement) println("------ORF " + orf_id + " not within bounds for " + x)
      statement
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
        if (!isWithinBoundaries(orf_i, x)) config.alpha
        else Math.abs(computeDisplacement(orf_i + x, orf_i) - phi(orf_i + x, orf_j, w))
      }
      //compute displacement in the perspective or orf i
      val perspective_orf_j = {
        if (!isWithinBoundaries(orf_j, y)) config.alpha
        else Math.abs(computeDisplacement(orf_j + y, orf_j) - phi(orf_j + y, orf_i, w))
      }
      if (config.verbose) println("----Found the following values " + (perspective_orf_i, perspective_orf_j))
      List(perspective_orf_i, perspective_orf_j).max * config.dpenalty
    }

    def constructSyntenicVectors: (Int, Int) =>
      (Seq[Double], Seq[Double], Seq[Double],Seq[Double]) = (orf_a, orf_b) => {
      if (config.verbose) println("--Using following positions away for left SV: " + (0 to config.f - 1))
      //compute left syntenic vector
      val (sv_l,sv_l_prime) = (0 to (config.f - 1)).foldLeft((Seq[Double](), Seq[Double]())){case((sv,svi), p) => {
        val (x,y) = (-(config.f - p), -(config.f - p))
        (sv.:+(psi(orf_a, orf_b, x, y, true)), svi.:+(psi(orf_a, orf_b, x, -y, false)))
      }}
      if (config.verbose) println("--Using following positions away for right SV: " + (1 to config.f))
      //compute left syntenic vector
      val (sv_r, sv_r_prime) = (1 to config.f).foldLeft((Seq[Double](),Seq[Double]())){case((sv,svi), p) => {
        val(x,y) = (p, p)
        (sv.:+(psi(orf_a, orf_b, x, y, true)),svi.:+(psi(orf_a, orf_b, x, -y,false)))
      }}
      if(config.verbose) {
        println("--Syntenic vectors: " + (sv_l,sv_r))
        println("--Inverse syntenic vectors: " + (sv_l_prime, sv_r_prime))
      }
      (sv_l, sv_r, sv_l_prime, sv_r_prime)
    }

    println("Computing syntenic scores for each BRH" + timeStamp)
    val syntenic_anchor_pairs = hashmap_H.toList.flatMap { case (orf, brhs) => {

      /**
        * Method that iterate through a list of BRHs and retains the best syntenic score along with corresponding ORF
        *
        * @param _brhs
        * @param min_score
        * @param min_orf
        * @return
        */
      def computeSyntenicScore(_brhs: List[Int],
                               min_score: Option[Double], min_orf: Option[Int]): (Option[Double], Option[Int]) = {
        if (_brhs.isEmpty) (min_score, min_orf)
        else {
          if(config.verbose) println("Scoring BRH: " + (orf, _brhs.head))
          //compute left and right syntenic vectors
          val (sv_l, sv_r, sv_l_prime, sv_r_prime) = constructSyntenicVectors(orf, _brhs.head)
          //compute synteny scores
          val (sv_l_score, sv_r_score, sv_l_prime_score, sv_r_prime_score) =
            (sv_l.sum, sv_r.sum, sv_l_prime.sum, sv_r_prime.sum)
          //compute syntenic score
          val syntenic_score = List(sv_l_score, sv_r_score, sv_l_prime_score, sv_r_prime_score).min
          if(config.verbose) println("--Syntenic score is: " + syntenic_score)
          pw_syntenic_scores.println(
            Seq(orf.toString, _brhs.head.toString, sv_l_score, sv_r_score, sv_l_prime_score, sv_r_prime_score).mkString("\t"))
          //update accordingly if the syntenic score is lower than previously observed
          if (min_score == None || syntenic_score < min_score.get){
            computeSyntenicScore(_brhs.tail, Option(syntenic_score), Option(_brhs.head))
          }
          else computeSyntenicScore(_brhs.tail, min_score, min_orf)
        }
      }

      //retain syntenic anchor per genome
      val pairs = brhs.groupBy(x => hashmap_Y_prime(hashmap_Z_prime(x))).mapValues(local_brhs => {
        //fetch best syntenic score and the corresponding orf
        val (best_syntenic_score, best_orf) = computeSyntenicScore(local_brhs, None, None)
        //control when there are unique sequences and when the syntenic scores is above specified threshold
        if (best_syntenic_score == None) (-1, -1)
        else if (!(best_syntenic_score.get <= (config.f * config.syntenicFraction))) (-1, -1)
        else {
          if(config.verbose) println("Classified as a syntenic anchor: " + (orf, best_orf))
          (orf, best_orf.get)
        }
      })
      //convert to list and eliminate cases of unique orfs
      pairs.toList.map(_._2).filterNot(_._1 == -1)
    }
    }
    //dump syntenic pairs
    if (config.dump) {
      val pw = new PrintWriter(config.outputDir + "/syntenic_pairs.txt")
      syntenic_anchor_pairs.foreach(x => pw.println(Seq(x._1, x._2).mkString("\t")))
      pw.close
    }
    println("Found " + syntenic_anchor_pairs.size + " syntenic anchor pairs" + timeStamp)
    val pw_syntenic_anchors = new PrintWriter(config.outputDir + "/syntenic_anchors.txt")
    println("Constructing syntenic anchors" + timeStamp)
    //construct syntenic anchors from syntenic pairs
    syntenic_anchor_pairs.foldLeft(HashMap[Int, Set[Int]]()) { case (map_syntenic_anchors, (orf_a, orf_b)) => {
      val current_a = map_syntenic_anchors.getOrElse(orf_a, Set[Int]())
      val current_b = map_syntenic_anchors.getOrElse(orf_b, Set[Int]())
      val update1 = (map_syntenic_anchors + (orf_a -> current_a.+(orf_b)))
      update1 + (orf_b -> current_b.+(orf_a))
    }//output to file
    }.foreach { case (orf, syntenic_anchors) => {
      pw_syntenic_anchors.println(Seq(orf, syntenic_anchors.mkString(",")).mkString("\t"))
    }
    }
    pw_syntenic_anchors.close
    pw_syntenic_scores.close
    println("Successfully completed!" + timeStamp)
  }
}
