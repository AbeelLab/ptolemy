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
      "Usage: java -jar ptolemy.jar [module]\n\n" +
        "CANONICAL QUIVER CONSTRUCTION\n" +
        "extract                      Run Ptolemy's extract module.\n"+
        "syntenic-anchors             Compute syntenic anchors from pairwise ORF alignments.\n"+
        "canonical-quiver             Construct canonical quiver.\n\n" +
        "STRUCTURAL VARIANT CALLING\n" +
        "variant-calling              Identify structural variants as a population using maximally-labelled paths\n\n"+
        "OPTIONAL\n" +
        "run-msa                      Induce MSA across syntenic anchors.\n"
      )

    if (args.length == 0) {
      println(help)
    } else {
      args(0) match {
        case "extract"             => build_db.Extract.main(args.drop(1))
        case "syntenic-anchors"    => syntenic_anchors.SyntenicAnchors.main(args.drop(1))
        case "canonical-quiver"    => canonical_quiver.ConstuctCanonicalQuiver.main(args.drop(1))
        case "run-msa"             => msa.RunMSA.main(args.drop(1))
        case "variant-calling"     => variant_calling.HighLevelVariantCaller.main(args.drop(1))
        case "parents"             => parental_perspective.CharacterizeWithParents.main(args.drop(1))
        case _                     => println(help)
      }
    }
  }

}
