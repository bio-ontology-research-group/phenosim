@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')
@Grab(group='org.codehaus.gpars', module='gpars', version='1.1.0')

import java.net.*
import org.openrdf.model.vocabulary.*
import slib.sglib.io.loader.*
import slib.sml.sm.core.metrics.ic.utils.*
import slib.sml.sm.core.utils.*
import slib.sglib.io.loader.bio.obo.*
import org.openrdf.model.URI
import slib.graph.algo.extraction.rvf.instances.*
import slib.sglib.algo.graph.utils.*
import slib.utils.impl.Timer
import slib.graph.algo.extraction.utils.*
import slib.graph.model.graph.*
import slib.graph.model.repo.*
import slib.graph.model.impl.graph.memory.*
import slib.sml.sm.core.engine.*
import slib.graph.io.conf.*
import slib.graph.model.impl.graph.elements.*
import slib.graph.algo.extraction.rvf.instances.impl.*
import slib.graph.model.impl.repo.*
import slib.graph.io.util.*
import slib.graph.io.loader.*
import groovyx.gpars.GParsPool

System.setProperty("jdk.xml.entityExpansionLimit", "0");
System.setProperty("jdk.xml.totalEntitySizeLimit", "0");

def factory = URIFactoryMemory.getSingleton()
def annotationsPath = "data/genes_to_phenotype.txt";
def resSimPath = "data/sim_gene_gene.txt";


def getOntology = {

  URI graph_uri = factory.getURI("http://purl.obolibrary.org/obo/")
  G graph = new GraphMemory(graph_uri)

  // Load OBO file to graph "go.obo"
  GDataConf goConf = new GDataConf(GFormat.RDF_XML, "data/hp.owl")
  GraphLoaderGeneric.populate(goConf, graph)

  // Add virtual root for 3 subontologies__________________________________
  URI virtualRoot = factory.getURI("http://purl.obolibrary.org/obo/virtualRoot")
  graph.addV(virtualRoot)

  new File(annotationsPath).splitEachLine('\t') { items ->
    if (items[0].startsWith("#")) return;
    String disId = items[0].replaceAll(":", "_");
    URI idURI = factory.getURI("http://purl.obolibrary.org/obo/" + disId);
    String pheno = items[3].replaceAll(":", "_");
    URI phenoURI = factory.getURI("http://purl.obolibrary.org/obo/" + pheno);
    Edge e = new Edge(idURI, RDF.TYPE, phenoURI);
    graph.addE(e);
  }

  GAction rooting = new GAction(GActionType.REROOTING)
  rooting.addParameter("root_uri", virtualRoot.stringValue())
  GraphActionExecutor.applyAction(factory, rooting, graph)
  return graph
}

def getURIfromName = { name ->
  // def id = name.split('\\:')
  return factory.getURI("http://purl.obolibrary.org/obo/$name")
}

graph = getOntology()
phenotypes = graph.getV();

def getMapping = {
  def map = [:].withDefault {new LinkedHashSet()}
  new File("data/9606.protein.aliases.v10.5.txt").splitEachLine('\t') { items ->
    if (items[0].startsWith("#")) return;
    def stId = items[0]
    def geneId = items[1]
    map[geneId].add(stId)
  }
  return map
}

def getGenes = {
  def mapping = getMapping()
  def genes = [:].withDefault {new LinkedHashSet()}
  new File(annotationsPath).splitEachLine('\t') { items ->
    if (items[0].startsWith("#")) return;
    def pheno = items[3].replaceAll(":", "_")
    def phenoURI = getURIfromName(pheno)
    if (!(phenoURI in phenotypes)) return;
    mapping[items[0]].each { geneId ->
      genes[geneId].add(phenoURI)
    }
  }
  return genes
}

def genes = getGenes()
def geneList = genes.keySet() as String[]
def n = geneList.size()
println("Number of genes: " + n)

def sim_id = 0 //this.args[0].toInteger()

SM_Engine engine = new SM_Engine(graph)

// BMA+Resnik, BMA+Schlicker2006, BMA+Lin1998, BMA+Jiang+Conrath1997,
// DAG-GIC, DAG-NTO, DAG-UI

String[] flags = [
  // SMConstants.FLAG_SIM_GROUPWISE_AVERAGE,
  // SMConstants.FLAG_SIM_GROUPWISE_AVERAGE_NORMALIZED_GOSIM,
  SMConstants.FLAG_SIM_GROUPWISE_BMA,
  SMConstants.FLAG_SIM_GROUPWISE_BMM,
  SMConstants.FLAG_SIM_GROUPWISE_MAX,
  SMConstants.FLAG_SIM_GROUPWISE_MIN,
  SMConstants.FLAG_SIM_GROUPWISE_MAX_NORMALIZED_GOSIM
]

// List<String> pairFlags = new ArrayList<String>(SMConstants.PAIRWISE_MEASURE_FLAGS);
String[] pairFlags = [
  SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995,
  SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_SCHLICKER_2006,
  SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998,
  SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_JIANG_CONRATH_1997_NORM
]

ICconf icConf = new IC_Conf_Corpus("ResnikIC", SMConstants.FLAG_IC_ANNOT_RESNIK_1995_NORMALIZED);
String flagGroupwise = flags[sim_id.intdiv(pairFlags.size())];
String flagPairwise = pairFlags[sim_id % pairFlags.size()];
SMconf smConfGroupwise = new SMconf(flagGroupwise);
SMconf smConfPairwise = new SMconf(flagPairwise);
smConfPairwise.setICconf(icConf);

// // Schlicker indirect
// ICconf prob = new IC_Conf_Topo(SMConstants.FLAG_ICI_PROB_OCCURENCE_PROPAGATED);
// smConfPairwise.addParam("ic_prob", prob);

def result = new Double[n][n]
def index = new Integer[n * n]
for (int i = 0; i < index.size(); i++) {
  index[i] = i
}

def c = 0

GParsPool.withPool {
  index.eachParallel { val ->
    def i = val.toInteger()
    def x = i.intdiv(n)
    def y = i % n
    if (x <= y) {
      result[x][y] = engine.compare(
            smConfGroupwise,
            smConfPairwise,
            genes[geneList[x]],
            genes[geneList[y]])
      result[y][x] = result[x][y]
      if (c % 100000 == 0)
        println c
      c++
    }
  }
}

def out = new PrintWriter(new BufferedWriter(
  new FileWriter(resSimPath)))
out.print(geneList[0]);
for (int i = 1; i < n; i++) {
  out.print("\t" + geneList[i]);
}
out.println();
for (int i = 0; i < n; i++) {
  out.print(result[i][0]);
  for (int j = 1; j < n; j++) {
      out.print("\t" + result[i][j]);
  }
  out.println();
}
out.flush()
out.close()
