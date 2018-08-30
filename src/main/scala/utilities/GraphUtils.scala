package utilities

import scala.annotation.tailrec
import scala.collection.immutable.HashMap

/**
  * Author: Alex N. Salazar
  * Created on 30-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GraphUtils {

  /**
    * Tail-recursive method to obtain all connected components in a given graph.
    *
    * @param graph          Graph as Map[Int, Set[Int]
    * @param global_visited All visited nodes
    * @param ccs            All connected components
    * @return List[Set[Int]
    */
  @tailrec def getConnectedComponents(graph: HashMap[Int, Set[Int]],
                             global_visited: Set[Int],
                             ccs: List[Set[Int]]): List[Set[Int]] = {
    /**
      * Tail-recursive method for constructing the complete graph given a set of ORF IDs
      *
      * @param queue
      * @param acc_ccs
      * @return
      */
    @tailrec def _getConnectedComponents(queue: Seq[Int],
                                locally_visited: Set[Int],
                                acc_ccs: Set[Int]): (Set[Int], Set[Int]) = {
      if (queue.isEmpty) (acc_ccs, locally_visited)
      else {
        //get brhs of current orf that have not been already visited
        val new_brhs = graph.getOrElse(queue.head, Set[Int]()).filterNot(locally_visited(_))
        //remove current from queue, add brhs to visited, and add current to acc and move to next
        _getConnectedComponents(queue.tail ++ new_brhs.toSeq,
          new_brhs.foldLeft(locally_visited)((b, a) => b + (a)),
          acc_ccs + queue.head)
      }
    }

    //no more nodes to process
    if (graph.isEmpty) ccs
    else {
      //get set of starting nodes
      val starting_nodes = graph.head._2 + graph.head._1
      val filtered_nodes = starting_nodes.filter(global_visited(_))
      //get size of starting nodes
      val s_size = starting_nodes.size
      //get size of starting nodes already been visited
      val f_size = filtered_nodes.size
      //sanity check: either all nodes have been processed, or none have. If not, failed to identify true complete
      // connected component
      assert(f_size == 0 || s_size - f_size == 0, "Starting sub-graph contains mixtures of nodes that been " +
        "unprocessed/processed:" + (starting_nodes, filtered_nodes))
      //if all nodes have been observed, move on
      if (f_size == s_size) getConnectedComponents(graph.tail, global_visited, ccs)
      //get connected component using these set of nodes
      else {
        //compute the complete graph for current orf
        val (local_cc, updated_visited) =
          _getConnectedComponents(starting_nodes.toSeq, global_visited ++ starting_nodes, Set())
        //remove all instances from the above that are still in the graph, append to group of complete graphs
        getConnectedComponents(graph.tail, updated_visited, ccs.:+(local_cc))
      }
    }
  }

}
