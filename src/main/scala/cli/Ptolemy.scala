package cli

/**
  * Author: Alex N. Salazar
  * Created on 13-12-2017
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object Ptolemy {

  def main(args: Array[String]): Unit = {
    val help = (
      "Usage: java -jar Pandora.jar [module]\n\n" +
        "extract                      Run Ptolemy's extract module.\n"+
        "syntenic-anchors             Compute syntenic anchors from pairwise ORF alignments.\n"+
        "construct-hlgg               Construct HLGG.\n"
      )

    if (args.length == 0) {
      println(help)
    } else {
      args(0) match {
        case "extract"             => extract_module.Extract.main(args.drop(1))
        case "syntenic-anchors"    => cluster_module.SyntenicAnchors.main(args.drop(1))
        case "construct-hlgg"      => construct.HLGG.main(args.drop(1))
        case _                     => println(help)
      }
    }
  }

}
