/*
 * This work is licensed under the Creative Commons
 * Attribution-NonCommercial-NoDerivs 3.0 Unported License.
 * To view a copy of this license, visit
 * http://creativecommons.org/licenses/by-nc-nd/3.0/
 * or send a letter to Creative Commons, 444 Castro Street,
 * Suite 900, Mountain View, California, 94041, USA.
 *
 * A copy of the license is included in LICENSE.txt
 * See the License for the specific language governing permissions
 * and limitations under the License.
 *
 * Copyright 2005-2016 Thomas Abeel
 */
package atk

import java.io.{File, PrintWriter}
import java.time.LocalDateTime

class NixWriter(f: String,config:AnyRef=null,noheader:Boolean=false) extends PrintWriter(f) with Tool{

  def this(f:File,config:AnyRef)={
    this(f.toString,config)
  }

  if(!noheader)
    println(generatorInfo(config))

  override def println(str: String) {
    print(str + "\n")
  }

  override def close(){
    println("## This analysis finished " + LocalDateTime.now())
    super.close

  }


}