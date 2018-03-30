package utilities

import java.io.{File, PrintWriter}

import quiver.{Context, Graph, LEdge, LNode, empty}
import utilities.FileHandling.openFileWithIterator
import utilities.GFAutils.ConstructGFA

import scala.annotation.tailrec
import scala.collection.immutable.{HashMap, HashSet}

/**
  * Author: Alex N. Salazar
  * Created on 27-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
class HLGG(gfa: File) extends ConstructGFA {

  /** Type alias for high level genome graph as a inductive graph */
  type HLGG = Graph[Int, Set[String], Set[String]]
  type HLGC = Context[Int, Set[String], Set[String]]
  val empty_HLGG = empty[Int, Set[String], Set[String]]


  /** Object storing HLGG */
  val (graph, sequence_set, genome_set) = contructHLGG
  val reverse_sequence_set = sequence_set.foldLeft(HashMap[Int, Set[String]]()){case (map, (seq, nodes)) => {
    nodes.foldLeft(map)((local_map, node) => {
      val get = local_map.getOrElse(node, Set[String]()) + seq
      local_map + (node -> get)
    })
  }}


  /**
    * Function to determine whether line is a path line
    *
    * @return Boolean
    */
  def isPathLine: String => Boolean = line => line.startsWith("P")

  def isGenomeLine: String => Boolean = line => line.startsWith("G")

  /**
    * Function to parse path line in GFA and return sequence identifier and path
    *
    * @return (String, Seq[Int])
    */
  private def getPath: String => (String, IndexedSeq[Int]) = line => {
    val split = line.split("\t")
    (split(1), split(2).split("""\W+""").map(_.toInt).toIndexedSeq)
  }
  /**
    * Function to parse path line in GFA and return sequence identifier and path
    *
    * @return (String, Seq[Int])
    */
  private def getGenomes: String => (String, Set[String]) = line => {
    val split = line.split("\t")
    (split(1), split(2).split(",").toSet)
  }


  /**
    * Function to update node in HLGG
    *
    * @return HLGG
    */
  def update4Node: (Int, Set[String], HLGG) => HLGG = (node, label, hlgg) => {
    if (!hlgg.contains(node)) hlgg.addNode(LNode(node, label))
    else {
      val current = hlgg.context(node)
      hlgg.updateNode(LNode(node, current.label ++ label))
    }
  }

  def update4Edge: (Int, Int, Set[String], HLGG) => HLGG = (n1, n2, labels, hlgg) => {
    val current = hlgg.context(n1).outEdges.find(_.to == n2)
    if (current == None) hlgg.addEdge(LEdge[Int, Set[String]](n1, n2, labels))
    else {
      hlgg.updateEdge(LEdge(n1, n2, current.get.label ++ labels))
    }
  }

  /**
    * Function to add a pair of nodes to HLGG
    *
    * @return HLGG
    */
  def insertNodes: (Int, Int, Set[String], HLGG) => HLGG = (n1, n2, label, hlgg) =>
    update4Node(n2, label, update4Node(n1, label, hlgg))

  /**
    * Function to add an edge to HLGG
    *
    * @return HLGG
    */
  def insertEdge: (Int, Int, Set[String], HLGG) => HLGG = (n1, n2, labels, hlgg) => update4Edge(n1, n2, labels, hlgg)


  /**
    * Function to add a pair of nodes and an their edge to HLGG
    *
    * @return
    */
  def addNodesWithEdges: (Int, Int, Set[String], HLGG) => HLGG = (n1, n2, labels, hlgg) => {
    insertEdge(n1, n2, labels, (insertNodes(n1, n2, labels, hlgg)))
  }


  /**
    * Function to construct HLGG, vertex and reverse vertex set given a GFA file. Will construct from path lines only.
    *
    * @return (HLGG, HahMap(Genome, Seq[Nodes]), HashMap(Node,
    */
  private def contructHLGG = {
    openFileWithIterator(gfa)
      .foldLeft((empty_HLGG, HashMap[String, IndexedSeq[Int]](), HashMap[String, String]())) {
        case ((hlgg, vs, gs), line) => {
      if (!isPathLine(line) && !isGenomeLine(line)) (hlgg, vs, gs)
      else if(!isPathLine(line) && isGenomeLine(line)){
        val (genome , seqs) = getGenomes(line)
        (hlgg, vs,
        seqs.foldLeft(gs)((local_gs, seq) => local_gs + (seq -> genome)))
      }
      else {
        val (seq_id, path) = getPath(line)
        path.size match {
          case 1 => (update4Node(path.head, Set(seq_id), hlgg), vs + (seq_id -> path), gs)
          //iterate through path as pair of nodes and add nodes and inferred edge to HLGG; also add hashmaps
          case _ => (path.sliding(2).foldLeft(hlgg)((lh, nodes) => {
              addNodesWithEdges(nodes(0), nodes(1), Set(seq_id), lh)
            }), vs + (seq_id -> path), gs)
        }
      }
    }
    }
  }

  def traverseExclusive(nodes: HashSet[Int]): HLGC => Seq[Int] = { context: HLGC =>
    context.outEdges.map(_.to).filter(!nodes(_))
  }

  def traverseReferenceFirst(genome: String): HLGC => Seq[Int] = { context: HLGC =>
    val outs = context.outs.map(x => graph.context(x._2)).filter(_.label == genome)
    if (!outs.isEmpty) outs.map(_.vertex) else context.outEdges.map(_.to)
  }

  def preOrderDFS(implicit hlgg: HLGG, unmarkedSet: Set[Int], vertextList: List[Int]): List[Int] = {
    vertextList match {
      case List() => Nil
      case head :: tail => {
        head :: preOrderDFS(hlgg, unmarkedSet - head,
          (hlgg.successors(head).toList.sorted ++ tail).intersect(unmarkedSet.toList))
      }
    }
  }

  /**
    * Method to traverse through graph using DFS and return post-order traversal
    *
    * @param unmarkedSet Set of all nodes in graph
    * @param vertexList  List containing single node to visit
    * @return List[Int] corresponding to post-order DFS traversal
    */
   def postOrderDFS(subgraph: HLGG)(unmarkedSet: Set[Int], vertexList: List[Int]): List[Int] = {
    vertexList match {
      case List() => Nil
      case head :: tail => {
        //Recursively visit successors that haven't been visited yet.
        //  Note: up to the point where there are no more children, will return empty list (see above). If so,
        //  appends the current node to the empty list. This corresponds to the finishing node in the current
        //  traversal
        val visited = postOrderDFS(subgraph)(unmarkedSet - head,
          graph.successors(head).sorted.toList.intersect(unmarkedSet.toList)).:+(head)
        //Since we are getting the post order, we add the current finishing node in front of the later finishing
        //  nodes acquired in the recursive call.
        visited ::: postOrderDFS(subgraph)(unmarkedSet -- visited, tail diff visited)
      }
    }
  }

  /**
    * Method to traverse through HLGG for each node given and return the maximal set of nodes obtained
    * through the traversals
    *
    * @param nodes
    * @param visited
    * @param traversals
    * @return
    */
  @tailrec final def collectDFStraversals(subgraph: HLGG)(nodes: List[Int], visited: HashSet[Int],
                                           traversals: List[List[Int]]): List[List[Int]] = {
    nodes match {
      case Nil => traversals
      case _ => {
        //traverse through graph using dfs on the head node; remove nodes in the traversal that were already visited
        val traversal = subgraph.dfs(nodes.take(1)).toList.filter(!visited(_))
        collectDFStraversals(subgraph)(nodes diff traversal, visited ++ traversal, traversal :: traversals)
      }
    }
  }

  /**
    * Method to obtain strong connected components using Kosraju's algorithm.
    *
    * @param subgraph Optionally provide list of nodes consisting
    * @return
    */
  def kosarajusSCC(subgraph: HLGG = empty_HLGG): List[List[Int]] = {
    val rpo = subgraph match {
      //use hlgg in class if subgraph not provided
      case x if (x.isEmpty) => postOrderDFS(graph)(graph.nodes.toSet, graph.nodes.toList.sorted)
      case _ => postOrderDFS(subgraph)(subgraph.nodes.toSet, subgraph.nodes.toList.sorted)
    }
    collectDFStraversals({
      if (subgraph.isEmpty) graph else subgraph
    })(rpo.reverse, HashSet(), List())
  }



}
