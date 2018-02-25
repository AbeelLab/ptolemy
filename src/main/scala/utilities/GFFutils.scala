package utilities

/**
  * Author: Alex N. Salazar
  * Created on 5-12-2017
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

trait GFFutils {

  /**Case class for GFFLine*/
  case class GFFLine(chrm: String, sample: String, feature: String,
                     start: Int, end: Int, orientation: String, name: Option[String], id: String) {
    //get parent attribute
    lazy val parent = if (feature == "gene") None else Option(id.replace("cds", "gene"))
    //get original attribute
    lazy val attribute =
      "ID=" + id +
        (if (parent != None) (";Parent=" + parent.get) else "") +
        (if (parent == None) ";Name=" + name.get else "")
    //override toString to return original GFFLine
    override def toString = {
      Seq(chrm, sample, feature, start, end, ".", orientation, ".", attribute).mkString("\t")
    }
  }

  /**
    * Function to parse gff line and return GFFLine object
    * @return GFFLine object
    */
  def toGFFLine: String => GFFLine = line => {
    val tmp = line.split("\t")
    new GFFLine(tmp(0), tmp(1), tmp(2), tmp(3).toInt, tmp(4).toInt, tmp(6), getGeneName(tmp(8)), getGeneID(tmp(8)))
  }

    /** Harcoded string in the GFF file to identify unique gene ID. Used in the method, 'getGeneID' */
    private val gene_id_string = "ID="

    /**
      * Method to parse out the gene name in the attributes column of a gff file
      */
    private def getGeneName(attribute: String): Option[String] = {
      val annotation = attribute.split(";")
      val locus_tag = annotation.find(_.contains("locus_tag="))
      val gene_name = annotation.find(_.contains("Name="))
      if (locus_tag == None && gene_name == None) None
      else if (locus_tag == None && gene_name != None) Option(gene_name.get.substring(5))
      else {
        if (locus_tag.get.substring(10) == gene_name.get.substring(5)) Option(locus_tag.get.substring(10))
        else Option(locus_tag.get.substring(10) + "_" + gene_name.get.substring(5))
      }
    }

    /**
      * Method to extract the local gene id
      */
    private def getGeneID(x: String): String = {
      val tmp = x.split(";").find(_.contains(gene_id_string)).get
      tmp.substring(tmp.indexOf(gene_id_string) + 3)
    }

}