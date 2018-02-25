package utilities

/**
  * Author: Alex N. Salazar
  * Created on 16-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
trait GFAutils {

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
    * @return
    */
  def addGenomeCountField: Option[Int] => String = genomes =>
    "\tFC:i:"+ {if(genomes == None) 1 else genomes.get}

  def addColourField: String => String = colour => "\tID:Z" + colour
}
