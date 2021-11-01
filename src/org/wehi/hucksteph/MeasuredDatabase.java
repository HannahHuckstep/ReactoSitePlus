package org.wehi.hucksteph;


import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.biopax.paxtools.model.level2.protein;
import org.neo4j.cypher.internal.frontend.v2_3.ast.functions.Str;
import org.neo4j.cypher.internal.frontend.v3_4.phases.Do;
import org.neo4j.graphalgo.GraphAlgoFactory;
import org.neo4j.graphalgo.PathFinder;
import org.neo4j.graphalgo.WeightedPath;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.graphdb.traversal.*;
import org.neo4j.unsafe.impl.batchimport.input.InputException;

import java.io.*;
import java.rmi.server.UID;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.Math.round;

public class MeasuredDatabase extends EmbeddedNeo4jDatabase{


    private final String UID_PATTERN = "[OPQ][0-9][A-Z0-9]{3}[0-9](\\-[0-9*]{1,2})?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(\\-[0-9*]{1,2})?";
    private final String humanUniProt = "https://www.uniprot.org/uniprot/?query=organism:9606&format=fasta&include=yes";
    private final String mouseUniProt = "https://www.uniprot.org/uniprot/?query=organism:10090&format=fasta&include=yes";


    public MeasuredDatabase(File databaseDir, File outputFile) {
        super(databaseDir, outputFile);
    }

    public MeasuredDatabase(File databaseDir) {
        super(databaseDir);
    }

    /**
     * @param node
     * @param graphDb
     * @return
     */
    private Traverser getDownstream(Node node, GraphDatabaseService graphDb) {

        TraversalDescription td = graphDb.traversalDescription().breadthFirst()
                .evaluator(new Evaluator() {
                    @Override
                    public Evaluation evaluate(Path path) {
                        if(path.endNode().getProperty(PropertyType.TYPE.toString())
                                .equals("SmallMolecule")){
                            return Evaluation.INCLUDE_AND_PRUNE;
                        }
                        return Evaluation.INCLUDE_AND_CONTINUE;
                    }
                })
                .relationships(RelTypes.PHOSPHORYLATION, Direction.INCOMING)
                .relationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING)
                .relationships(RelTypes.INPUT, Direction.OUTGOING)
                .relationships(RelTypes.OUTPUT, Direction.OUTGOING)
                .relationships(RelTypes.CONTROLS, Direction.OUTGOING)
                .relationships(RelTypes.CATALYSIS, Direction.OUTGOING)
                .relationships(RelTypes.COMPONENT, Direction.OUTGOING)
                .relationships(RelTypes.SMALL_MOL_EDGE, Direction.OUTGOING)
                .relationships(RelationshipType.withName("ACTIVATION"), Direction.OUTGOING)
                .relationships(RelationshipType.withName("INHIBITION"), Direction.OUTGOING);
        return td.traverse(node);
    }

    /**
     * @param node
     * @param graphDb
     * @return
     */
    private Traverser getUpstream(Node node, GraphDatabaseService graphDb) {
        TraversalDescription td = graphDb.traversalDescription().breadthFirst()
                .evaluator(new Evaluator() {
                    @Override
                    public Evaluation evaluate(Path path) {
                        if(path.startNode().getProperty(PropertyType.TYPE.toString())
                                .equals("SmallMolecule")){
                            return Evaluation.INCLUDE_AND_PRUNE;
                        }
                        return Evaluation.INCLUDE_AND_CONTINUE;
                    }
                })
                .relationships(RelTypes.INPUT, Direction.INCOMING)
                .relationships(RelTypes.OUTPUT, Direction.INCOMING)
                .relationships(RelTypes.CONTROLS, Direction.INCOMING)
                .relationships(RelTypes.PHOSPHORYLATION, Direction.INCOMING)
                .relationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING)
                .relationships(RelTypes.COMPONENT, Direction.INCOMING)
                .relationships(RelTypes.SMALL_MOL_EDGE, Direction.INCOMING)
                .relationships(RelationshipType.withName("ACTIVATION"), Direction.INCOMING)
                .relationships(RelationshipType.withName("INHIBITION"), Direction.INCOMING);
        return td.traverse(node);
    }

    /**
     * Overloaded getRXNs mehtod that takes in a file as an input
     * @param pathToProts The path to a file containing a list of uniprot IDs or DB_ID's to be traversed
     * @param direction The direction of the traversal ("Upstream" or "Downstream")
     * @throws IOException
     */
    public void traversal(File pathToProts, String direction, String experiment) throws IOException {
        // get Proteins from File
        BufferedReader BR = null;
        String line = "";
        String csvSplitBy = ",";
        BR = new BufferedReader(new FileReader(pathToProts));
        Set<String> UIDList = new HashSet<>();
        while ((line = BR.readLine()) != null) {
            UIDList.add(line);
        }

        for (String uid : UIDList) {
            traversal(uid, direction, experiment);
        }
    }

    /**
     * Finds everything up or downstream of a node of interest
     * @param uid the starting UID or DB_ID to br traversed
     * @param direction the direction of the traversal ("Upstream" or "Downstream")
     * @param experiment the name of the experiment the statistics will be reported on
     * @throws IOException
     */
    public void traversal(String uid, String direction, String experiment) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        File OUTPUT_PATH = outputFile;

        direction = direction.toLowerCase();

        FileWriter fstream3 = new FileWriter(OUTPUT_PATH + "/R_input_traversal_lengths_"+uid+"_"+direction+".tsv");
        BufferedWriter out3 = new BufferedWriter(fstream3);
        out3.write("NodeID\texperiment\tUID\tdispName\tUniProtName\tlocation\tmoleculesDS\tpf\tpfMapped\tphos\tphosMapped\tcplx\tcplxMapped\tkins\tkinsMapped\tUids\tUidsMapped\tbchmRxns\n");

        try (Transaction tx = graphDb.beginTx()) {
            Pattern p = Pattern.compile(UID_PATTERN);
            Matcher m = p.matcher(uid);
            if (m.find()) { //////////////////////////// if its a uniprotid
                String theGroup = m.group(0);
                Node node = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()),
                        PropertyType.UNIPROT_ID.toString(), theGroup);

                // get all of its uid nodes
                // if it doesnt exist continue
                if (node == null) {
                    System.out.println(theGroup + " Does not exist in the databse");
                } else {
                    FileWriter fstream = new FileWriter(OUTPUT_PATH + "/"+uid+"_"+direction+"_traversal_cytoscape.tsv");
                    BufferedWriter out = new BufferedWriter(fstream);
                    //out.write("ID_"+node.getId()+"\t" + uid + "_"+direction+"\n");
                    out.write("ID_"+node.getId()+"\t");

                    FileWriter fstream1 = new FileWriter(OUTPUT_PATH + "/TraversalReport_"+direction +"_"+uid+".tsv");
                    BufferedWriter out1 = new BufferedWriter(fstream1);

                    // store traversers to sort later
                    HashMap<Node, Traverser> traverserMap = new HashMap<>();
                    Iterable<Relationship> relationships = node.getRelationships(RelTypes.ID_BELONGS_TO);
                    for (Relationship relationship : relationships) {
                        Node prot = relationship.getEndNode();
                        if (direction.equalsIgnoreCase("downstream")) {

                            Traverser nodeTraverser = getDownstream(prot, graphDb);
                            traverserMap.put(prot, nodeTraverser);

                        } else if (direction == "upstream") {

                            Traverser nodeTraverser = getUpstream(prot, graphDb);
                            traverserMap.put(prot, nodeTraverser);

                        } else {
                            throw new IllegalArgumentException("direction must equal 'upstream' or 'downstream'");
                        }
                    }


                    // get lengths of each traverser
                    HashMap<Node, Integer> traverserLen = new HashMap<>();
                    HashSet<Integer> uniqueLens = new HashSet<>();
                    for (Node key: traverserMap.keySet()) {
                        Traverser paths = traverserMap.get(key);
                        int len = getTraverserLen(paths);
                        traverserLen.put(key, len);
                        uniqueLens.add(len);
                    }

                    // now you have 2 maps Node:Traverser, Node:TraverserLen
                    // write to stream

                    // first get unique set of all things ds of all initial pf's
                    // and make hashmap of all streams to their nodes -> Node:[dsNode1, dsNode2]
                    HashMap<Node, HashSet<Node>> nodeMap = new HashMap<>();
                    HashSet<Node> uniqueAll = new HashSet<>();
                    for (Node key: traverserMap.keySet()) {
                        Traverser paths = traverserMap.get(key);
                        ResourceIterable<Node> alldsnodes = paths.nodes();
                        HashSet<Node> alldsnodesSet = new HashSet<>();
                        for (Node n: alldsnodes) {
                            uniqueAll.add(n);
                            alldsnodesSet.add(n);
                        }
                        nodeMap.put(key, alldsnodesSet);
                    }

                    for (Node key: nodeMap.keySet()) {
                        out.write( key.getId() + "_stream\t");
                    }
                    out.write("\n");

                    for (Node n: uniqueAll) {
                        String s = n.getId() + "\t";
                        for (Node key: nodeMap.keySet()) {
                            if (nodeMap.get(key).contains(n)){
                                s = s.concat(String.valueOf(key.getId()));
                            }else{
                                s = s.concat("");
                            }
                            s = s.concat("\t");
                        }
                        s = s.concat("\n");
                        out.write(s);
                    }


                    // print largest, delete largest from list
                    Collection<Integer> values = traverserLen.values();
                    List<Integer> lens = new ArrayList<>(values);
                    Integer count =  0;

                    Integer max = Collections.max(lens);
                    if(values.size() >1){
                        for (int i = 0; i < uniqueLens.size(); i++) {
                            count++;
                            for (Node key: traverserLen.keySet()) {
                                if (traverserLen.get(key).equals(max)){
                                    traversalReport(count, traverserMap.get(key), graphDb,uid,key,out1, experiment, out3);
                                }
                            }
                            lens.removeAll(Collections.singleton(max));
                            if(!lens.isEmpty()){
                                max = Collections.max(lens);
                            }
                        }
                    }else{
                        for (Node key: traverserLen.keySet()) {
                            traversalReport(count, traverserMap.get(key), graphDb,uid,key,out1, experiment, out3);
                        }

                    }

                    out.close();
                    out1.close();
                }
            } else { ///////////////// its a nodeID
                try{
                    Long.valueOf(uid);
                }catch (InputException e){
                    System.out.println("Input must be a UniProt id (ex. P04637), a node id (ex. 12), or a file of UniProt ids");
                    e.printStackTrace();
                    System.exit(1);
                }

                Node prot = graphDb.getNodeById(Long.valueOf(uid));
                // get all of its uid nodes
                // if it doesnt exist continue

                if (direction == "downstream") {
                    FileWriter fstream = new FileWriter(OUTPUT_PATH + "/"+uid+"_Downstream.tsv");
                    BufferedWriter out = new BufferedWriter(fstream);
                    out.write("ID\t" + uid + "_downstream\n");

                    FileWriter fstream1 = new FileWriter(OUTPUT_PATH + "/TraversalReport_downstream_"+uid+".tsv");
                    BufferedWriter out1 = new BufferedWriter(fstream1);

                    Traverser nodeTraverser = getDownstream(prot, graphDb);
                    traversalReport(0, nodeTraverser, graphDb, uid, prot, out1, experiment, out3);

                    ResourceIterable<Node> nodes = nodeTraverser.nodes();
                    for (Node n: nodes) {
                        out.write(n.getId() + "\t" + prot.getId() + "\n");
                    }

                    out.close();
                    out1.close();

                } else if (direction == "upstream") {

                    FileWriter fstream2 =  new FileWriter(OUTPUT_PATH + "/"+uid+"_Upstream.tsv");
                    BufferedWriter in  = new BufferedWriter(fstream2);
                    in.write("ID\t" + uid + "_upstream\n");

                    FileWriter fstream1 = new FileWriter(OUTPUT_PATH + "/TraversalReport_upstream_"+uid+".tsv");
                    BufferedWriter out1 = new BufferedWriter(fstream1);

                    Traverser nodeTraverser = getUpstream(prot, graphDb);
                    traversalReport(0,nodeTraverser, graphDb, uid, prot, out1, experiment,out3);

                    ResourceIterable<Node> nodes = nodeTraverser.nodes();
                    for (Node n: nodes) {
                        in.write(n.getId() + "\t" + prot.getId() + "\n");
                    }

                    in.close();
                    out1.close();

                } else {
                    throw new IllegalArgumentException("direction must equal 'upstream' or 'downstream'");

                }
            }

            //}
            tx.success();
            out3.close();
        }

        graphDb.shutdown();
    }

    private Integer getTraverserLen(Traverser nodeTraverser){

        Integer depth = 0;
        for (Path path: nodeTraverser) {
            depth++;
        }
        return depth;
    }

    private void traversalReport(Integer order,
                                 Traverser nodeTraverser,
                                 GraphDatabaseService graphDb,
                                 String uid,
                                 Node prot,
                                 BufferedWriter reportStream,
                                 String experiment,
                                 BufferedWriter travLenStream) {

        File databaseDir = getDatabaseDir();

        try (Transaction tx = graphDb.beginTx()) {

            // get all labels and find one that matches the experiment
            Boolean experimentLabel = true;
            String supportScoreString = "";
            HashSet<String> experiments = new HashSet<>();
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            for (String property: allPropertyKeys){
                if(property.replace("SUPPORT_SCORE_", "").equalsIgnoreCase(experiment)){ // if experiment name is in properties
                    experimentLabel = false; // dont throw error
                    supportScoreString = property;
                }else if (property.contains("SUPPORT_SCORE_")){ // gather experiment names in db
                    String support_score_ = property.replace("SUPPORT_SCORE_", "");
                    experiments.add(support_score_);
                }
            }
            // throw exception if experiment name given is not in the database
            if(experimentLabel){
                throw new InputException(experiment + " is not currently in this database: " + databaseDir +
                        "\nExperiments in this database are: " + experiments);

            }


            int depth = 0;
            Integer phosCounter = 0;
            Integer phosMappedCounter = 0;
            Integer proteinCounter = 0;
            Integer proteinMappedCounter = 0;
            Integer cplxCounter = 0;
            Integer cplxMappedCounter = 0;
            Integer uidCounter = 0;
            Integer uidMappedCounter = 0;
            Integer transcriptionFactorCounter = 0;
            Integer cellSurfaceReceptorCounter = 0;
            Integer kinaseCounter = 0;
            Integer kinaseMappedCounter = 0;
            Integer rxnCounter = 0;
            HashSet<String> pathwaysTraversed = new HashSet<>();



            for (Path nodePath : nodeTraverser) {
                depth++;
                //phos
                if (nodePath.endNode().hasLabel(Label.label(LabelTypes.PHOSPHORYLATION.toString()))){
                    phosCounter ++;
                    Iterable<Relationship> relationships = nodePath.endNode().getRelationships(RelTypes.PHOSPHORYLATION);
                    for (Relationship relationship : relationships) {
                        Node protein = relationship.getEndNode();
                        if(protein.hasProperty(supportScoreString)){
                            phosMappedCounter ++;
                        }
                    }
                }
                //protein
                else if (nodePath.endNode().hasLabel(Label.label("Protein"))){
                    proteinCounter ++;
                    if(nodePath.endNode().hasProperty(supportScoreString)){
                        proteinMappedCounter++;
                    }
                    if(nodePath.endNode().hasProperty(PropertyType.KINASE.toString())){
                        kinaseCounter++;
                        if(nodePath.endNode().hasProperty(supportScoreString)){
                            kinaseMappedCounter++;
                        }
                    }
                    if(nodePath.endNode().hasProperty(PropertyType.TRANSCRIPTION_FACTOR.toString())){
                        transcriptionFactorCounter++;
                    }
                    if(nodePath.endNode().hasProperty(PropertyType.CELL_SURFACE_RECEPTOR.toString())){
                        cellSurfaceReceptorCounter++;
                    }
                }
                //complex
                else if (nodePath.endNode().hasLabel(Label.label("Complex"))){
                    cplxCounter ++;
                    if(nodePath.endNode().hasProperty(supportScoreString)){
                        cplxMappedCounter++;
                    }
                }
                //uid
                else if(nodePath.endNode().hasLabel(Label.label(LabelTypes.UNIPROT_ID.toString()))){
                    uidCounter++;
                    Iterable<Relationship> relationships = nodePath.endNode().getRelationships(RelTypes.ID_BELONGS_TO);
                    Boolean uidMappedTo = false;
                    for (Relationship relationship : relationships) { // dont want to count all proteins with meausrements attached to this uid
                        Node protein = relationship.getEndNode();
                        if(protein.hasProperty(supportScoreString)){
                            uidMappedTo = true; // so if this UID has at least 1 measured protein count just it\
                        }
                    }
                    if(uidMappedTo){
                        uidMappedCounter++;
                    }
                }
                //rxn
                else if (nodePath.endNode().hasLabel(Label.label("BiochemicalReaction"))){
                    rxnCounter++;
                }
                // pathways traversed
                Iterable<Relationship> relationships = nodePath.endNode().getRelationships(RelTypes.PATHWAY_COMPONENT);
                for (Relationship relationship : relationships) {
                    String property = relationship.getStartNode().getProperty(PropertyType.DISPLAY_NAME.toString()).toString();
                    pathwaysTraversed.add(property);
                }


            }
            String geneName = "";
            if(prot.hasProperty(PropertyType.UNIPROT_NAME.toString())){
                geneName = prot.getProperty(PropertyType.UNIPROT_NAME.toString()).toString();
            }
            travLenStream.write(prot.getId() +"\t"+
                    experiment +"\t"+
                    uid +"\t"+
                    prot.getProperty(PropertyType.DISPLAY_NAME.toString()) +"\t"+
                    geneName +"\t"+
                    prot.getProperty(PropertyType.LOCATION.toString())+"\t"+
                    depth +"\t"+
                    proteinCounter+"\t"+
                    proteinMappedCounter+"\t"+
                    phosCounter+"\t"+
                    phosMappedCounter+"\t"+
                    cplxCounter+"\t"+
                    cplxMappedCounter+"\t"+
                    kinaseCounter+"\t"+
                    kinaseMappedCounter+"\t"+
                    uidCounter+"\t"+
                    uidMappedCounter+"\t"+
                    rxnCounter+"\n");

            try {
                reportStream.write("Order: "+ order +"\nAll things Downstream of " + uid +
                        ", gene name:"+geneName+
                        ", display name: "+ prot.getProperty(PropertyType.DISPLAY_NAME.toString()) +
                        " (location: "+prot.getProperty(PropertyType.LOCATION.toString())+
                        ")(node id: "+prot.getId()+
                        ") in experiment " + experiment);
                reportStream.write("\nThe number of nodes found were: "+depth);
                reportStream.write("\nThe number of phosphorylation nodes: " + phosCounter);
                reportStream.write("\nThe number of phosphorylation nodes with mapped data: " + phosMappedCounter +" ("+proportion(phosCounter, phosMappedCounter)+"%)");
                reportStream.write("\nThe number of proteoform nodes: " + proteinCounter);
                reportStream.write("\nThe number of proteoform nodes with mapped data: " + proteinMappedCounter +"("+proportion(proteinCounter, proteinMappedCounter)+"%)");
                reportStream.write("\nThe number of complex nodes: " + cplxCounter);
                reportStream.write("\nThe number of complex nodes with mapped data: " + cplxMappedCounter +"("+proportion(cplxCounter, cplxMappedCounter)+"%)");
                reportStream.write("\nThe number of UniProt id nodes: " + uidCounter);
                reportStream.write("\nThe number of UniProt id nodes with mapped data: " + uidMappedCounter +"("+proportion(uidCounter, uidMappedCounter)+"%)");
                reportStream.write("\nThe number of kinase nodes: " + kinaseCounter);
                reportStream.write("\nThe number of kinase nodes with mapped data: " + kinaseCounter +"("+proportion(kinaseCounter, kinaseMappedCounter)+"%)");
                reportStream.write("\nThe number of biochemical reaction nodes: " + rxnCounter);
                reportStream.write("\nThe number of transcription factor nodes: " + transcriptionFactorCounter);
                reportStream.write("\nThe number of cell surface receptor nodes: " + cellSurfaceReceptorCounter);
                reportStream.write("\nThe number of pathways traversed: " + pathwaysTraversed.size());
                /*
                reportStream.write("\n\nThe pathways traversed: ");
                for (String path: pathwaysTraversed) {
                    reportStream.write("\n" + path);
                }
                 */
                reportStream.write("\n\n");

            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            }

            tx.success();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


    }

    private Double proportion(Integer num, Integer numMapped){
        Double total = Double.valueOf(num);
        Double mapped = Double.valueOf(numMapped);

        Double prop = mapped/total;
        prop = (Math.round(prop*100.00))/1.00;


        return prop;
    }

    private Traverser getNeighbourhoodBFS(Node node, Integer depth, GraphDatabaseService graphDb) {

        TraversalDescription td = graphDb.traversalDescription().breadthFirst()
                .evaluator(new Evaluator() {
                    @Override
                    public Evaluation evaluate(Path path) {
                        if(path.endNode().getProperty(PropertyType.TYPE.toString())
                                .equals("SmallMolecule")){
                            return Evaluation.INCLUDE_AND_PRUNE;
                        }
                        return Evaluation.INCLUDE_AND_CONTINUE;
                    }
                })
                .evaluator(Evaluators.toDepth(depth))
                .relationships(RelTypes.OUTPUT, Direction.BOTH)
                .relationships(RelTypes.INPUT, Direction.BOTH)
                .relationships(RelTypes.CONTROLS, Direction.BOTH)
                .relationships(RelTypes.CATALYSIS, Direction.BOTH)
                .relationships(RelTypes.SMALL_MOL_EDGE, Direction.BOTH)
                .relationships(RelationshipType.withName("ACTIVATION"), Direction.BOTH)
                .relationships(RelationshipType.withName("INHIBITION"), Direction.BOTH);

        return td.traverse(node);
    }

    /**
     * actually runs the traversal and returns the list of protein nodes in neighbourhood
     * @param node
     * @param depth
     * @param graphDb
     * @return
     */
    HashSet<Node> neighbourTraversal(Node node, Integer depth, GraphDatabaseService graphDb){
        Integer currdepth = 0;

        String output = "starting at: " + node.getId() + "\n";
        Traverser nodeTraverser = getNeighbourhoodBFS(node, depth, graphDb);
        HashSet<Node> nbhd = new HashSet<>();
        for (Path nodePath : nodeTraverser) {
            output += "at depth" + nodePath.length() + "=>" +nodePath.endNode().getLabels()+" " + nodePath.endNode().getId() + "\n";
            currdepth++;
            //if(nodePath.endNode().hasLabel(Label.label("Protein"))){
            nbhd.add(nodePath.endNode());
            //}
        }
        output += "Number of nodes found: " + currdepth + "\n";
        //System.out.println(output);

        return nbhd;
    }

    /**
     * Performs neighbourhood analysis using the biomial distribution
     * @param depth
     * @throws IOException
     */
    public void binomialNeighbourhood(Integer depth) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);


        try(Transaction tx = graphDb.beginTx()){
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            List<String> scores = new ArrayList<>();
            for (String property: allPropertyKeys) {
                if (property.startsWith("SUPPORT_SCORE_")){
                    scores.add(property);
                }
            }
            FileWriter fstream3 = new FileWriter(outputFile + "/BinomialNeighbourhoodsReport.tsv");
            BufferedWriter out3 = new BufferedWriter(fstream3);

            FileWriter fstream4 = new FileWriter(outputFile + "/experiments.tsv");
            BufferedWriter out4 = new BufferedWriter(fstream4);
            out4.write("Experiments\tdepth");
            // for each experiment type, do a neighbourhood analysis
            for (String scoreString: scores) {

                if(scoreString.contains("SUPPORT_SCORE_") | scoreString.contains("ABUNDANCE_SCORE_") ){
                    scoreString = scoreString.replaceAll("SUPPORT_SCORE_", "");
                    scoreString = scoreString.replaceAll("ABUNDANCE_SCORE_", "");
                }
                out3.write("\n\nExperiment: " + scoreString);
                out4.write("\n"+scoreString + "\t" + depth);
                System.out.println("\nCurently analysing: " + scoreString);
                String supportScoreStr = "SUPPORT_SCORE_" + scoreString;
                String abundScoreStr = "ABUNDANCE_SCORE_" + scoreString;
                String scoredByStr = "SCORED_BY_" + scoreString;

                // write new files
                FileWriter fstream = new FileWriter(outputFile + "/BinomialNeighbourhoods_"+scoreString + ".tsv");
                BufferedWriter out = new BufferedWriter(fstream);
                out.write("ID\tDispName\tUID\tPval\tBonfCorrected\tAvgSScore\tAvgAScore\tNumInNbhd\tNumMeasuredNbhd\tNumIntegrated\tNumUIDsInNBHD\tdbMeasured\tprobablility\n");

                FileWriter fstream2 = new FileWriter(outputFile + "/SigNeighbourhoods_"+scoreString+".tsv");
                BufferedWriter out2 = new BufferedWriter(fstream2);
                out2.write("ID\tDispName\tUID\tPval\tBonfCorrected\tAvgSScore\tAvgAScore\tNumInNbhd\tNumMeasuredNbhd\tNumIntegrated\tNumUIDsInNBHD\tdbMeasured\tprobablility\n");


                // get the number of measured proteins and complexes
                int dbMeasured = 0;
                double numInDb = 0;
                ResourceIterator<Node> allPEs = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
                while (allPEs.hasNext()){
                    Node pe = allPEs.next();
                    if(pe.hasLabel(Label.label("Protein")) | pe.hasLabel(Label.label("Complex"))){
                        if (pe.hasProperty(supportScoreStr)){
                            dbMeasured++;
                        }
                        numInDb++;
                    }
                }


                Integer count = 0;
                Integer numBonfSig = 0;
                Integer numPvalSig = 0;
                // for each protein and complex generate the neighbourhood and get it's pvalue
                ResourceIterator<Node> PEs = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
                //ResourceIterator<Node> PEs = graphDb.findNodes(Label.label("Protein"));
                for (ResourceIterator<Node> it = PEs; it.hasNext(); ) {
                    count ++;
                    ArrayList<Double> nbhdSScores = new ArrayList<>();
                    ArrayList<Double> nbhdAScores = new ArrayList<>();

                    if(count % 10000 == 0 ){
                        System.out.print("\rProgress: "+ (count/numInDb)*100);
                    }

                    Node pe = it.next();
                    if(pe.hasLabel(Label.label("Protein")) | pe.hasLabel(Label.label("Complex"))){

                        HashSet<Node> nbhd = neighbourTraversal(pe, depth, graphDb);

                        HashSet<String> nbhdUids = new HashSet<>();

                        // count the total neighbourhood size
                        // count number of measured things
                        double nbhdSize = 0;
                        int nbhdMeasured = 0;
                        int numIntegrated = 0;
                        for(Node node:nbhd){
                            if (node.hasProperty(supportScoreStr)){
                                if(Double.parseDouble(node.getProperty(supportScoreStr).toString()) > 0){
                                    nbhdMeasured ++;
                                }
                                nbhdSScores.add(Double.parseDouble(node.getProperty(supportScoreStr).toString()));
                                nbhdAScores.add(Double.parseDouble(node.getProperty(abundScoreStr).toString()));
                                nbhdSize ++;

                                if(node.hasProperty(PropertyType.INTEGRATED.toString())){
                                    numIntegrated++;
                                }

                                // count num unique uids in nbhd
                                Iterable<Relationship> uidRels = node.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                                if (getLength(uidRels) == 1){
                                    for(Relationship uidRel: uidRels){
                                        nbhdUids.add(uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                                    }
                                }else{
                                    String uids = "";
                                    for(Relationship uidRel: uidRels){
                                        uids = uids + uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString().toString());
                                    }
                                    nbhdUids.add(uids);
                                }
                            }
                        }



                        double probability = nbhdSize/numInDb;

                        Double pval = 0.0;

                        if(nbhdMeasured != 0){
                            BinomialTest bn = new BinomialTest();
                            pval = bn.binomialTest(dbMeasured, nbhdMeasured, probability, AlternativeHypothesis.GREATER_THAN);
                        }

                        Double bonfCorrected= pval * dbMeasured;


                        ////////////// gathering outputs
                        String uid = "";
                        if(pe.hasProperty(scoredByStr)){
                            uid = pe.getProperty(scoredByStr).toString();
                        } else if (pe.hasProperty(PropertyType.UNIPROT_ID.toString())){
                            uid = pe.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                        }else{
                            Iterable<Relationship> uidRels = pe.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                            if(getLength(uidRels) == 1){
                                for(Relationship uidRel: uidRels){
                                    uid = uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                                }
                            }else{
                                for(Relationship uidRel: uidRels){
                                    uid = uid + uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString() + ",";
                                }
                            }
                        }

                        Double ssAvg = 0.0;
                        Double asAvg = 0.0;

                        for(Integer i = 0; i < nbhdSScores.size(); i++){
                            if(nbhdSScores.get(i) >= 0){
                                ssAvg += nbhdSScores.get(i);
                                asAvg += nbhdAScores.get(i);
                            }
                        }

                        ssAvg = ssAvg/nbhdSScores.size();
                        asAvg = asAvg/nbhdSScores.size();


                    /*
                    if(pval == 0){
                        System.out.println(pe.getProperty(PropertyType.DISPLAY_NAME.toString())
                                + ", Num in DB: " + dbMeasured
                                + ", NBHD measured: " + nbhdMeasured
                                +  ", prob: " + probability
                                + " = " + nbhdSize +"/"+ numInDb
                        );
                    }
                     */

                        out.write( pe.getId() + "\t" +
                                pe.getProperty(PropertyType.DISPLAY_NAME.toString()) + "\t" +
                                uid + "\t" +
                                pval + "\t" +
                                bonfCorrected + "\t" +
                                ssAvg+"\t"+
                                asAvg + "\t"+
                                nbhdSize + "\t" +
                                nbhdMeasured  + "\t" +
                                numIntegrated  + "\t" +
                                nbhdUids.size() + "\t" +
                                dbMeasured + "\t" +
                                probability + "\n"
                        );

                        if(pval < 0.05){
                            numPvalSig ++;
                        }

                        if(bonfCorrected < 0.05 & nbhdSize > 0.0 & ssAvg > 0.0){
                            numBonfSig ++;
                            out2.write( pe.getId() + "\t" +
                                    pe.getProperty(PropertyType.DISPLAY_NAME.toString()) + "\t" +
                                    uid + "\t" +
                                    pval + "\t" +
                                    bonfCorrected + "\t" +
                                    ssAvg+"\t"+
                                    asAvg + "\t"+
                                    nbhdSize + "\t" +
                                    nbhdMeasured  + "\t" +
                                    numIntegrated  + "\t" +
                                    nbhdUids.size() + "\t" +
                                    dbMeasured + "\t" +
                                    probability + "\n"
                            );
                        }
                    }
                }

                out3.write("\nThe number of neighbourhoods tested: " + numInDb);
                out3.write("\nThe number of neighbourhoods with a p-value < 0.05: " + numPvalSig);
                out3.write("\nThe number of neighbourhoods with a Bonferroni corrected p-value < 0.05: " + numBonfSig);

                tx.success();
                out.close();
                out2.close();

            }
            out3.close();
            out4.close();

        }
        graphDb.shutdown();
    }

    /**
     * Uses dijkstra to find the shortest path between 2 nodes, can be Uniprot ID's or node ids and path can be upstream or downstream
     * @param startString
     * @param endString
     * @throws IOException
     */
    public void shortestPath(String startString, String endString, String weightType) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        FileWriter fstream = new FileWriter(outputFile + "/ShortestPathReport_"+startString+"_to_"+endString+".tsv");
        BufferedWriter out = new BufferedWriter(fstream);
        out.write("Shortest Paths Report: ");


        try (Transaction tx = graphDb.beginTx()) {
            // get all labels and find one that matches the experiment
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            for (String property: allPropertyKeys){
                // get experiment name for weights
                String experimentStr = "";
                String weightString = "";
                String supportScoreStr = "";
                if(property.contains("SUPPORT_SCORE_")){
                    // get each experiments weight and support  property
                    supportScoreStr = property;
                    experimentStr = property.replaceAll("SUPPORT_SCORE_" , "");

                    if(weightType.equalsIgnoreCase("Abundance") | weightType.equalsIgnoreCase("a")){
                        weightString = "WEIGHT_ABUNDANCE_" + experimentStr;
                    }else if(weightType.equalsIgnoreCase("Support") | weightType.equalsIgnoreCase("s")){
                        weightString = "WEIGHT_SUPPORT_" + experimentStr;
                    }

                    Node start = null;
                    Node end = null;
                    Pattern p = Pattern.compile(UID_PATTERN);
                    Matcher m = p.matcher(startString);
                    if (m.find()) { // if its a uniprotid
                        String theGroup = m.group(0);
                        start = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()),
                                PropertyType.UNIPROT_ID.toString(), theGroup);
                        if (start == null) {
                            throw new NotFoundException(startString+" does not exist in the database " + databaseDir);
                        }
                    }else{ // its a node id
                        start = graphDb.getNodeById(Long.valueOf(startString));
                        if (start == null) {
                            throw new NotFoundException(startString+" does not exist in the database " + databaseDir);
                        }
                    }

                    Matcher m1 = p.matcher(endString);
                    if (m1.find()) { // if its a uniprotid
                        String theGroup = m1.group(0);
                        end = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()),
                                PropertyType.UNIPROT_ID.toString(), theGroup);
                        if (end == null) {
                            throw new NotFoundException(endString+" does not exist in the database " + databaseDir);
                        }
                    }else{// its a node id
                        end = graphDb.getNodeById(Long.valueOf(endString));
                        if (end == null) {
                            throw new NotFoundException(endString+" does not exist in the database " + databaseDir);
                        }
                    }

                    PathFinder<WeightedPath> dijkstraDS = GraphAlgoFactory.dijkstra(
                            PathExpanders.forTypesAndDirections(
                                    RelTypes.INPUT, Direction.OUTGOING,
                                    RelTypes.OUTPUT, Direction.OUTGOING,
                                    RelTypes.CONTROLS, Direction.OUTGOING,
                                    RelTypes.CATALYSIS, Direction.OUTGOING,
                                    RelTypes.ID_BELONGS_TO, Direction.OUTGOING,
                                    RelTypes.ID_BELONGS_TO, Direction.INCOMING,
                                    RelTypes.COMPONENT, Direction.OUTGOING,
                                    RelationshipType.withName("ACTIVATION"), Direction.OUTGOING,
                                    RelationshipType.withName("INHIBITION"), Direction.OUTGOING),
                            weightString,
                            50
                    );

                    PathFinder<WeightedPath> dijkstraUS = GraphAlgoFactory.dijkstra(
                            PathExpanders.forTypesAndDirections(
                                    RelTypes.INPUT, Direction.INCOMING,
                                    RelTypes.OUTPUT, Direction.INCOMING,
                                    RelTypes.CONTROLS, Direction.INCOMING,
                                    RelTypes.CATALYSIS, Direction.INCOMING,
                                    RelTypes.ID_BELONGS_TO, Direction.OUTGOING,
                                    RelTypes.ID_BELONGS_TO, Direction.INCOMING,
                                    RelTypes.COMPONENT, Direction.INCOMING,
                                    RelationshipType.withName("ACTIVATION"), Direction.INCOMING,
                                    RelationshipType.withName("INHIBITION"), Direction.INCOMING),
                            weightString,
                            50
                    );

                    // find and report the shortest path stats for each experiment mapped

                    Iterable<WeightedPath> allPathsDS = dijkstraDS.findAllPaths(start, end);
                    Iterator<WeightedPath> iteratorDS = allPathsDS.iterator();
                    Integer lengthDS = getLength(iteratorDS);
                    Iterable<WeightedPath> allPathsUS = dijkstraUS.findAllPaths(start, end);
                    Iterator<WeightedPath> iteratorUS = allPathsUS.iterator();
                    Integer lengthUS = getLength(iteratorUS);
                    if(lengthDS == 0 & lengthUS == 0){
                        throw new Exception("No path between nodes "+startString+" and "+endString +" either upstream or downstream");
                    } if (lengthDS > 0){

                        FileWriter fstream1 = new FileWriter(outputFile + "/"+startString+"_to_"+endString+"_downstream_cytoscape.tsv");
                        BufferedWriter out1 = new BufferedWriter(fstream1);
                        out1.write("nodeID\t"+startString+"_to_"+endString+"\n");


                        for (Path path: allPathsDS){ // DOWNSTREAM
                            //get details
                            Integer numMolecules = 0;
                            Integer numProteins = 0;
                            Integer numMeasuredProteins = 0;
                            Integer numPhosdProteins = 0;
                            Integer numMeasuredPhosdProteins = 0;
                            Integer numComplexes = 0;
                            Integer numMeasuredComplexes = 0;
                            Integer numRelationships = 0;
                            Iterable<Node> nodes = path.nodes();
                            Iterator<Node> iterator = nodes.iterator();
                            while (iterator.hasNext()){
                                Node node = iterator.next();
                                // store last 3 node types in array
                                // if pattern is UPU
                                    // break


                                if(node.hasLabel(Label.label("Protein"))){
                                    numMolecules++;
                                    numProteins++;
                                    out1.write(node.getId()+"\t"+startString+"_to_"+endString+"\n");
                                    if(node.hasProperty(supportScoreStr)){
                                        numMeasuredProteins++;
                                    }
                                    if(node.hasRelationship(RelTypes.PHOSPHORYLATION)){
                                        numPhosdProteins++;
                                        if (node.hasProperty(supportScoreStr)){
                                            numMeasuredPhosdProteins++;
                                        }
                                    }
                                }else if(node.hasLabel(Label.label("Complex"))){
                                    numMolecules++;
                                    numComplexes++;
                                    out1.write(node.getId()+"\t"+startString+"_to_"+endString+"\n");
                                    if (node.hasProperty(supportScoreStr)){
                                        numMeasuredComplexes++;
                                    }
                                }else if(node.hasLabel(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()))){
                                    numMolecules++;
                                    out1.write(node.getId()+"\t"+startString+"_to_"+endString+"\n");
                                }else if(node.hasLabel(Label.label(LabelTypes.INTERACTION.toString()))){
                                    numRelationships++;
                                    out1.write(node.getId()+"\t"+startString+"_to_"+endString+"\n");
                                }
                            }

                            // print 'em
                            out.write("\n\nPath between " + startString +" and " + endString + " in experiment " + supportScoreStr.replace("SUPPORT_SCORE_", ""));
                            out.write("\n" + startString +" (" +start.getId()+ ") is upstream of " + endString +" ("+end.getId()+")" );
                            out.write("\nLength of path: " + path.length());
                            out.write("\nNumber of molecules: " + numMolecules);
                            out.write("\nNumber of reactions: " + numRelationships);
                            out.write("\nNumber of proteins: " + numProteins);
                            out.write("\nNumber of measured proteins: " + numMeasuredProteins);
                            out.write("\nNumber of phosphorylated proteins: " + numPhosdProteins);
                            out.write("\nNumber of measured phosphorylated proteins: " + numMeasuredPhosdProteins);
                            out.write("\nNumber of complexes: " + numComplexes);
                            out.write("\nNumber of measured complexes: " + numMeasuredComplexes);
                            out.write("\nPath with id's:\n" + path.toString());
                            out.write("\nPath with names:\n");
                            Iterator<PropertyContainer> iterator1 = path.iterator();
                            while (iterator1.hasNext()){
                                PropertyContainer next = iterator1.next();
                                if(next instanceof Node){
                                    if(next.hasProperty(PropertyType.DISPLAY_NAME.toString())){
                                        out.write(next.getProperty(PropertyType.DISPLAY_NAME.toString()).toString() + "->");
                                    }else{
                                        long id = ((Node) next).getId();
                                        out.write("("+String.valueOf(id) + ")->");
                                    }
                                }else if(next instanceof Relationship){
                                    RelationshipType type = ((Relationship) next).getType();
                                    if(next.hasProperty(weightString)){
                                        String weight = next.getProperty(weightString).toString();
                                        out.write("["+type+"]("+weight+")->");
                                    }else{
                                        out.write("["+type+"]->");
                                    }
                                }
                            }
                        }
                        out1.close();
                    }else if (lengthUS > 0){

                        FileWriter fstream1 = new FileWriter(outputFile + "/"+startString+"_to_"+endString+"_upstream_cytoscape.tsv");
                        BufferedWriter out1 = new BufferedWriter(fstream1);
                        out1.write("nodeID\t"+startString+"_to_"+endString+"\n");

                        for (Path path: allPathsUS){ // UPSTREAM
                            //get details
                            Integer numMolecules = 0;
                            Integer numProteins = 0;
                            Integer numMeasuredProteins = 0;
                            Integer numPhosdProteins = 0;
                            Integer numMeasuredPhosdProteins = 0;
                            Integer numComplexes = 0;
                            Integer numMeasuredComplexes = 0;
                            Integer numRelationships = 0;
                            Iterable<Node> nodes = path.nodes();
                            Iterator<Node> iterator = nodes.iterator();
                            while (iterator.hasNext()){
                                Node node = iterator.next();
                                if(node.hasLabel(Label.label("Protein"))){
                                    numMolecules++;
                                    numProteins++;
                                    out1.write(node.getId()+"\t"+startString+"_to_"+endString+"\n");
                                    if(node.hasProperty(supportScoreStr)){
                                        numMeasuredProteins++;
                                    }
                                    if(node.hasRelationship(RelTypes.PHOSPHORYLATION)){
                                        numPhosdProteins++;
                                        if (node.hasProperty(supportScoreStr)){
                                            numMeasuredPhosdProteins++;
                                        }
                                    }
                                }else if(node.hasLabel(Label.label("Complex"))){
                                    numMolecules++;
                                    numComplexes++;
                                    out1.write(node.getId()+"\t"+startString+"_to_"+endString+"\n");
                                    if (node.hasProperty(supportScoreStr)){
                                        numMeasuredComplexes++;
                                    }
                                }else if(node.hasLabel(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()))){
                                    numMolecules++;
                                    out1.write(node.getId()+"\t"+startString+"_to_"+endString+"\n");
                                }else if(node.hasLabel(Label.label(LabelTypes.INTERACTION.toString()))){
                                    numRelationships++;
                                    out1.write(node.getId()+"\t"+startString+"_to_"+endString+"\n");
                                }
                            }

                            // print 'em
                            out.write("\n\nPath between " + startString +" and " + endString + " in experiment " + supportScoreStr.replace("SUPPORT_SCORE_", ""));
                            out.write("\n"+ startString +" ("+start.getId()+ ") is downstream of " + endString +" ("+end.getId()+")");
                            out.write("\nLength of path: " + path.length());
                            out.write("\nNumber of molecules: " + numMolecules);
                            out.write("\nNumber of reactions: " + numRelationships);
                            out.write("\nNumber of proteins: " + numProteins);
                            out.write("\nNumber of measured proteins: " + numMeasuredProteins);
                            out.write("\nNumber of phosphorylated proteins: " + numPhosdProteins);
                            out.write("\nNumber of measured phosphorylated proteins: " + numMeasuredPhosdProteins);
                            out.write("\nNumber of complexes: " + numComplexes);
                            out.write("\nNumber of measured complexes: " + numMeasuredComplexes);
                            out.write("\nPath with id's:\n" + path.toString());
                            out.write("\nPath with names:\n");
                            Iterator<PropertyContainer> iterator1 = path.iterator();
                            while (iterator1.hasNext()){
                                PropertyContainer next = iterator1.next();
                                if(next instanceof Node){
                                    if(next.hasProperty(PropertyType.DISPLAY_NAME.toString())){
                                        out.write("<-"+next.getProperty(PropertyType.DISPLAY_NAME.toString()).toString());
                                    }else{
                                        long id = ((Node) next).getId();
                                        out.write("<-("+String.valueOf(id) + ")");
                                    }
                                }else if(next instanceof Relationship){
                                    RelationshipType type = ((Relationship) next).getType();
                                    if(next.hasProperty(weightString)){
                                        String weight = next.getProperty(weightString).toString();
                                        out.write("<-["+type+"]("+weight+")");
                                    }else{
                                        out.write("<-["+type+"]");
                                    }
                                }
                            }
                        }
                        out1.close();
                    }
                }
            }

            tx.success();
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        graphDb.shutdown();

    }

    public void minimalConnectionNetwork(String experiment) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        try (Transaction tx = graphDb.beginTx()) {
            // get all labels and find one that matches the experiment
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            HashSet<String> experimentStrs = new HashSet<>();
            Boolean measured = false;
            for (String property : allPropertyKeys) {
                // get experiment name for weights
                if(property.replace("SUPPORT_SCORE_", "").equalsIgnoreCase(experiment)){ // if experiment name is in properties
                    measured = true;
                }else if (property.contains("SUPPORT_SCORE_")){ // gather experiment names in db
                    String support_score_ = property.replace("SUPPORT_SCORE_", "");
                    experimentStrs.add(support_score_);
                }
            }
            if(!measured){
                throw new InputException(experiment + " is not currently in this database: " + databaseDir +
                        "\nExperiments in this database are: " + experimentStrs);
            }

            System.out.println("\nMCN FOR "+ experiment);
            FileWriter fstream = new FileWriter(outputFile + "/MinimalConnectionNetworkReport_"+experiment+".tsv");
            BufferedWriter out = new BufferedWriter(fstream);
            out.write("Minimal Connection Network (MCN) Report: ");
            out.write("\nExperiment: " + experiment);

            String supportScoreStr = "SUPPORT_SCORE_" + experiment;
            /*
            String weightString = "";
            if(weightType.equalsIgnoreCase("Abundance") | weightType.equalsIgnoreCase("a")){
                weightString = "WEIGHT_ABUNDANCE_" + experiment;
            }else if(weightType.equalsIgnoreCase("Support") | weightType.equalsIgnoreCase("s")){
                weightString = "WEIGHT_SUPPORT_" + experiment;
            }*/

            HashSet<Node> mapped = new HashSet<>();
            ResourceIterator<Node> nodes = graphDb.findNodes(Label.label("Protein"));
            while (nodes.hasNext()){
                Node node = nodes.next();
                if(node.hasProperty(supportScoreStr)){
                    mapped.add(node);
                }
            }
            System.out.println("Mapped Proteoforms: "+ mapped.size());

            PathExpanderBuilder builder = PathExpanderBuilder.allTypes(Direction.OUTGOING)
                    .remove(RelationshipType.withName(PropertyType.SMALL_MOL_EDGE.toString()))
                    .remove(RelTypes.PATHWAY_COMPONENT) // always incoming, but just to be safe
                    .remove(RelTypes.ID_BELONGS_TO); // always incoming, but just to be safe
            PathFinder<Path> finder = GraphAlgoFactory.shortestPath(builder.build(), 10);

            // make a hashset for all nodes and all relationships
            HashSet<Node> mcnNodes = new HashSet<>();
            HashSet<Relationship> mcnRels = new HashSet<>();

            Integer count = 0;
            int mappedSize = mapped.size();

            // find paths between all
            for (Node start: mapped){

                HashSet<Integer> uidsStart = new HashSet<Integer>();
                count ++;
                System.out.print("\rProgress: "+ (count/mappedSize)*100);

                Iterable<Relationship> uidRels = start.getRelationships(RelTypes.ID_BELONGS_TO);
                for(Relationship rel: uidRels){
                    uidsStart.add((int) rel.getStartNode().getId());
                }

                mcnNodes.add(start);
                for(Node end: mapped){

                    HashSet<Integer> uidsEnd = new HashSet<Integer>();
                    Iterable<Relationship> uidRels2 = end.getRelationships(RelTypes.ID_BELONGS_TO);
                    for(Relationship rel: uidRels2){
                        uidsEnd.add((int) rel.getStartNode().getId());
                    }
                    // if not the same node and nodes dont have same UIDs
                    if(start.getId() != end.getId() & Collections.disjoint(uidsStart, uidsEnd)){
                        //System.out.println("finding path between: "+ start.getId() + "and"+ end.getId());
                        //Iterable<WeightedPath> allPathsDS = dijkstraDS.findAllPaths(start, end);
                        //Iterator<WeightedPath> iterator = allPathsDS.iterator();
                        Iterable<Path> allPaths = finder.findAllPaths(start, end);
                        Iterator<Path> iterator = allPaths.iterator();
                        while (iterator.hasNext()){
                            //WeightedPath path = iterator.next();

                            Path path = iterator.next();
                            for (Node node : path.nodes()) {
                                mcnNodes.add(node);
                            }
                            for (Relationship relationship : path.relationships()) {
                                mcnRels.add(relationship);
                            }
                            //mcnNodes.addAll((Collection<? extends Node>) path.nodes());
                            //mcnRels.addAll((Collection<? extends Relationship>) path.relationships());
                        }
                    }
                }
            }

            FileWriter fstream1 = new FileWriter(outputFile + "/MinimalConnectionNetwork_"+experiment+"_cytoscape.tsv");
            BufferedWriter out1 = new BufferedWriter(fstream1);
            out1.write("NodeID\tMCN_"+experiment+"\n");

            Integer numProteins = 0;
            Integer numMappedProteins = 0;
            Integer numPhosdProteins = 0;
            Integer numMappedPhosdProteins = 0;
            Integer numComplexes = 0;
            Integer numMappedComplexes = 0;
            HashSet<String> numUIDs = new HashSet<>();
            Integer numRXNs = 0;
            Integer numKinases = 0;
            Integer numMappedKinases = 0;
            Integer numIntegrated = 0;

            for (Node node: mcnNodes){
                out1.write(node.getId() + "\tMCN_"+experiment+"\n");
                if(node.hasLabel(Label.label("Protein"))){
                    numProteins++;
                    if(node.hasProperty(supportScoreStr)){
                        numMappedProteins++;
                    }
                    if(node.hasProperty(PropertyType.INTEGRATED.toString())){
                        numIntegrated++;
                    }
                    if (node.hasRelationship(RelTypes.PHOSPHORYLATION)){
                        numPhosdProteins++;
                        if(node.hasProperty(supportScoreStr)){
                            numMappedPhosdProteins++;
                        }
                    }
                    if(node.hasProperty(PropertyType.KINASE.toString())) {
                        numKinases++;
                        if (node.hasProperty(supportScoreStr)) {
                            numMappedKinases++;
                        }
                    }
                    Iterable<Relationship> uidRels = node.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                    for(Relationship uidRel: uidRels){
                        Node uidNode = uidRel.getStartNode();
                        numUIDs.add(uidNode.getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                        out1.write(uidNode.getId() + "\tMCN_"+experiment+"\n");
                    }

                }else if(node.hasLabel(Label.label("Complex"))){
                    numComplexes++;
                    if(node.hasProperty(supportScoreStr)){
                        numMappedComplexes++;
                    }
                    if(node.hasProperty(PropertyType.INTEGRATED.toString())){
                        numIntegrated++;
                    }
                }else if (node.hasLabel(Label.label("BiochemicalReaction"))){
                    numRXNs++;
                }
            }

            for (String uid: numUIDs) {
                out1.write(uid + "\tMCN_"+experiment+"\n");
            }
            out1.close();

            out.write("\nThe number of proteoforms in the MCN: " + numProteins);
            out.write("\nThe number of mapped proteofroms in the MCN: " + numMappedProteins);
            out.write("\nThe number of phosphorylated proteoforms in the MCN: " + numPhosdProteins);
            out.write("\nThe number of mapped phosphorylated proteoforms in the MCN: " + numMappedPhosdProteins);
            out.write("\nThe number of complexes in the MCN: " + numComplexes);
            out.write("\nThe number of mapped complexes in the MCN: "+ numMappedComplexes);
            out.write("\nThe number of UniProt ids in the MCN: "+ numUIDs.size());
            out.write("\nThe number of mapped biochemical reactions in the MCN: "+numRXNs);
            out.write("\nThe number of kinases in the MCN: "+numKinases);
            out.write("\nThe number of mapped kinases in the MCN: "+numMappedKinases);
            out.write("\nThe number of integrated nodes in the MCN: "+numIntegrated);
            out.write("\n");

            out.close();
            tx.success();
        } catch (Exception e) {
            e.printStackTrace();
        }
        graphDb.shutdown();
        System.out.println();
    }

    /**
     * Takes a graph and removes any properties that have a score associated with the experiment name given
     */
    void resetScores( String experiment){
        File databaseDir = getDatabaseDir();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        try (Transaction tx = graphDb.beginTx()) {

            // get all labels and find one that matches the experiment
            Boolean experimentLabel = true;
            String supportScoreString = "";
            String abundanceScoreString = "";
            String mappedString = "";
            String weightSString = "";
            String weightAString = "";
            String scoredByString = "";
            HashSet<String> experiments = new HashSet<>();
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            for (String property: allPropertyKeys){
                if(property.replace("SUPPORT_SCORE_", "").equalsIgnoreCase(experiment)){ // if experiment name is in properties
                    experimentLabel = false; // dont throw error
                    supportScoreString = property;
                    abundanceScoreString = property.replaceAll("SUPPORT_SCORE", "ABUNDANCE_SCORE");
                    mappedString = property.replaceAll("SUPPORT_SCORE", "MAPPED");
                    weightAString = property.replaceAll("SUPPORT_SCORE", "WEIGHT_ABUNDANCE");
                    weightSString = property.replaceAll("SUPPORT_SCORE", "WEIGHT_SUPPORT");
                    scoredByString = property.replaceAll("SUPPORT_SCORE_", "SCORED_BY");
                }else if (property.contains("SUPPORT_SCORE_")){ // gather experiment names in db
                    String support_score_ = property.replace("SUPPORT_SCORE_", "");
                    experiments.add(support_score_);
                }
            }
            // throw exception if experiment name given is not in the database
            if(experimentLabel){
                throw new Exception(experiment + " is not currently in this database: " + databaseDir +
                        "\nExperiments in this database are: " + experiments);
            }

            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            for (ResourceIterator<Node> it = proteins; it.hasNext(); ) {
                Node protein = it.next();
                protein.removeProperty(supportScoreString);
                protein.removeProperty(abundanceScoreString);
                protein.removeProperty(scoredByString);
            }

            ResourceIterator<Node> complexs = graphDb.findNodes(Label.label("Complex"));
            for (ResourceIterator<Node> it = complexs; it.hasNext(); ) {
                Node complex = it.next();
                complex.removeProperty(supportScoreString);
                complex.removeProperty(abundanceScoreString);
                complex.removeProperty(scoredByString);
            }

            ResourceIterator<Node> uids = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            for (ResourceIterator<Node> it = uids; it.hasNext(); ) {
                Node uid = it.next();
                uid.removeProperty(mappedString);
            }

            for (Relationship allRelationship : graphDb.getAllRelationships()) {
                allRelationship.removeProperty(weightAString);
                allRelationship.removeProperty(weightSString);

            }

            for (Node allNode : graphDb.getAllNodes()) {
                if (allNode.hasProperty(supportScoreString)){
                    System.out.println(allNode.getAllProperties());
                }
            }


            tx.success();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Takes a graph and removes any properties that have a score associated with the experiment name given
     */
    void resetScores( GraphDatabaseService graphDb, String experiment){
        File databaseDir = getDatabaseDir();

        try (Transaction tx = graphDb.beginTx()) {

            // get all labels and find one that matches the experiment
            Boolean experimentLabel = true;
            String supportScoreString = "";
            String abundanceScoreString = "";
            String mappedString = "";
            String weightSString = "";
            String weightAString = "";
            String scoredByString = "";
            HashSet<String> experiments = new HashSet<>();
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            for (String property: allPropertyKeys){
                if(property.replace("SUPPORT_SCORE_", "").equalsIgnoreCase(experiment)){ // if experiment name is in properties
                    experimentLabel = false; // dont throw error
                    supportScoreString = property;
                    abundanceScoreString = property.replaceAll("SUPPORT_SCORE", "ABUNDANCE_SCORE");
                    mappedString = property.replaceAll("SUPPORT_SCORE", "MAPPED");
                    weightAString = property.replaceAll("SUPPORT_SCORE", "WEIGHT_ABUNDANCE");
                    weightSString = property.replaceAll("SUPPORT_SCORE", "WEIGHT_SUPPORT");
                    scoredByString = property.replaceAll("SUPPORT_SCORE_", "SCORED_BY");
                }else if (property.contains("SUPPORT_SCORE_")){ // gather experiment names in db
                    String support_score_ = property.replace("SUPPORT_SCORE_", "");
                    experiments.add(support_score_);
                }
            }
            // throw exception if experiment name given is not in the database
            if(experimentLabel){
                throw new Exception(experiment + " is not currently in this database: " + databaseDir +
                        "\nExperiments in this database are: " + experiments);
            }

            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            for (ResourceIterator<Node> it = proteins; it.hasNext(); ) {
                Node protein = it.next();
                protein.removeProperty(supportScoreString);
                protein.removeProperty(abundanceScoreString);
                protein.removeProperty(scoredByString);
            }

            ResourceIterator<Node> complexs = graphDb.findNodes(Label.label("Complex"));
            for (ResourceIterator<Node> it = complexs; it.hasNext(); ) {
                Node complex = it.next();
                complex.removeProperty(supportScoreString);
                complex.removeProperty(abundanceScoreString);
                complex.removeProperty(scoredByString);
            }

            ResourceIterator<Node> uids = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            for (ResourceIterator<Node> it = uids; it.hasNext(); ) {
                Node uid = it.next();
                uid.removeProperty(mappedString);
            }

            for (Relationship allRelationship : graphDb.getAllRelationships()) {
                allRelationship.removeProperty(weightAString);
                allRelationship.removeProperty(weightSString);

            }

            tx.success();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private Integer getLength(Iterable<Relationship> thing) {
        Integer count = 0;
        for (Relationship relationship : thing) {
            count++;
        }
        return count;
    }

    private Integer getLength(Iterator<WeightedPath> thing) {
        Integer count = 0;
        while(thing.hasNext()){
            Object object = thing.next();
            count ++;
        }
        return count;
    }

    void empiricalNullDistribution(File qPhosFile, Integer depth, Integer subsetSize, Integer repetitionNumber ) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        // reset scores
        try (Transaction tx = graphDb.beginTx()) {
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            for (ResourceIterator<Node> it = proteins; it.hasNext(); ) {
                Node protein = it.next();
                protein.removeProperty(PropertyType.SUPPORT_SCORE.toString());
                protein.removeProperty(PropertyType.ABUNDANCE_SCORE.toString());
                protein.removeProperty(PropertyType.SCORED_BY.toString());
            }

            ResourceIterator<Node> complexs = graphDb.findNodes(Label.label("Complex"));
            for (ResourceIterator<Node> it = complexs; it.hasNext(); ) {
                Node complex = it.next();
                complex.removeProperty(PropertyType.SUPPORT_SCORE.toString());
                complex.removeProperty(PropertyType.ABUNDANCE_SCORE.toString());
                complex.removeProperty(PropertyType.SCORED_BY.toString());
                complex.removeProperty(PropertyType.MAPPED.toString());
            }

            ResourceIterator<Node> uids = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            for (ResourceIterator<Node> it = uids; it.hasNext(); ) {
                Node uid = it.next();
                uid.removeProperty(PropertyType.MAPPED.toString());
            }
            tx.success();
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Reading file in and getting hashmaps for mapping

        // Read in data into Dict
        BufferedReader BR = new BufferedReader(new FileReader(qPhosFile));
        String line = "";

        HashMap<String,String> ID2uniquePep = new HashMap<>();
        HashMap<String, HashSet<String>> ID2LRP = new HashMap<>();
        HashMap<String, HashSet<Integer>> ID2Positions = new HashMap<>();
        HashMap<String, HashSet<String>> ID2TopMods = new HashMap<>();
        HashMap<String, Double> ID2aScore = new HashMap<>();
        HashMap<String, String> ID2RawPep = new HashMap<>();
        HashMap<String, String> ID2Start = new HashMap<>();
        HashMap<String, String> ID2End = new HashMap<>();
        HashSet<String> ids2Remove = new HashSet<>();
        // b/c the peptide has multi phosn's but less than that have recorded info so we dont know what mod in the peptide is the one with the measurements


        Integer leadingProtCol =0;
        Integer modPositionCol = 0;
        Integer aScoreCol = 0;
        Integer rawPepCol = 0;
        Integer conditionCol = 0;

        System.out.println("Reading in all qPhos Data");

        Integer countID = 0;
        while ((line = BR.readLine()) != null) {
            String[] columns = line.split("\t");
            Integer colLen = columns.length;
            if(countID == 0){ // if first line in file
                for (int i = 0; i < colLen; i++) {
                    if (columns[i].equalsIgnoreCase("UniProt accession")){
                        leadingProtCol = i;
                    }
                    if(columns[i].equals("Position")){
                        modPositionCol = i;
                    }
                    if(columns[i].equals("Log2Ratio")){
                        aScoreCol = i;
                    }
                    if(columns[i].equals("Raw peptide")){
                        rawPepCol = i;
                    }
                    if(columns[i].equals("Condition")){
                        conditionCol = i;
                    }
                }
                countID++;
            }else{ //for all other lines

                String UID_new = columns[leadingProtCol];
                String rawPep_new = columns[rawPepCol];
                String aScore_new = columns[aScoreCol];
                String position_new = columns[modPositionCol];
                String condition_new = columns[conditionCol];

                String keyStr = UID_new+rawPep_new+aScore_new+condition_new;

                /*
                System.out.println();
                System.out.println(keyStr);
                System.out.println("ID2uniquePep: " + ID2uniquePep);
                System.out.println("ID2LRP: " + ID2LRP);
                System.out.println("ID2Pos: " + ID2Positions);
                System.out.println("ID2RawPep: " + ID2RawPep);
                 */

                if(ID2uniquePep.containsKey(keyStr)){ // if you've seen this pep before
                    // get it's ID and add to it's dicts
                    String ID = ID2uniquePep.get(keyStr);
                    HashSet<Integer> positions = ID2Positions.get(String.valueOf(ID));
                    positions.add(Integer.valueOf(position_new));
                    ID2Positions.put(String.valueOf(ID), positions);

                }else{ // if this is a new entry (pep)
                    ID2uniquePep.put(keyStr, String.valueOf(countID));
                    // create new entries to all dicts
                    HashSet<Integer> positions = new HashSet<>();

                    if(ID2LRP.containsKey(UID_new)){
                        HashSet<String> temp = ID2LRP.get(UID_new);
                        temp.add(String.valueOf(countID));
                        ID2LRP.put(UID_new, temp);
                    }else{
                        HashSet<String> temp = new HashSet<>();
                        temp.add(String.valueOf(countID));
                        ID2LRP.put(UID_new, temp);
                    }

                    ID2aScore.put(String.valueOf(countID), Double.valueOf(aScore_new));
                    ID2RawPep.put(String.valueOf(countID), rawPep_new);
                    positions.add(Integer.valueOf(position_new));
                    ID2Positions.put(String.valueOf(countID), positions);

                    // incr ID #
                    countID ++;
                }
            }
        }



        // now you have:
        // ID2Positions -> ID: 4,7 (single ID to multi mod positions)
        // ID2RawPep -> ID: KDsGLYLKELIEPVLtCFNDADsRL

        for(String ID: ID2Positions.keySet()){

            String rawPep = ID2RawPep.get(ID);
            List<Integer> positions = new ArrayList<>(ID2Positions.get(ID));
            Collections.sort(positions);

            // so for each position, create a mod string ( S_4, T_7)
            Integer modCount = 0;
            for (int i = 0; i <rawPep.length() ; i++) {



                if(rawPep.charAt(i) == 's'){
                    if(modCount > (positions.size()-1)){
                        ids2Remove.add(ID);
                        continue;
                    }
                    Integer position = positions.get(modCount);
                    if(ID2TopMods.containsKey(ID)){
                        HashSet<String> temp = ID2TopMods.get(ID);
                        temp.add("S_"+position);
                        ID2TopMods.put(ID, temp);
                    }else{
                        HashSet<String> temp = new HashSet<>();
                        temp.add("S_" + position);
                        ID2TopMods.put(ID, temp);
                    }
                    modCount ++;

                } else if(rawPep.charAt(i) == 't'){
                    if(modCount > (positions.size()-1)){
                        ids2Remove.add(ID);
                        continue;
                    }
                    Integer position = positions.get(modCount);
                    if(ID2TopMods.containsKey(ID)){
                        HashSet<String> temp = ID2TopMods.get(ID);
                        temp.add("T_"+position);
                        ID2TopMods.put(ID, temp);
                    }else{
                        HashSet<String> temp = new HashSet<>();
                        temp.add("T_" + position);
                        ID2TopMods.put(ID, temp);
                    }
                    modCount ++;

                } else if(rawPep.charAt(i) == 'y'){
                    if(modCount > (positions.size()-1)){
                        ids2Remove.add(ID);
                        continue;
                    }
                    Integer position = positions.get(modCount);
                    if(ID2TopMods.containsKey(ID)){
                        HashSet<String> temp = ID2TopMods.get(ID);
                        temp.add("Y_"+position);
                        ID2TopMods.put(ID, temp);
                    }else{
                        HashSet<String> temp = new HashSet<>();
                        temp.add("Y_" + position);
                        ID2TopMods.put(ID, temp);
                    }
                    modCount ++;
                }
            }

            // so for each position, calculate a  start and end of pep
            // use first position
            Integer position = positions.get(0);
            Integer start = 0;
            Integer end = 0;

            for (int i = 0; i < rawPep.length() ; i++) {
                if(rawPep.charAt(i) == 's' | rawPep.charAt(i) == 't' | rawPep.charAt(i) == 'y'){
                    start = position - i;
                    end = start + rawPep.length()-1;
                    break;
                }
            }
            ID2Start.put(ID, String.valueOf(start));
            ID2End.put(ID, String.valueOf(end));
        }

        /*HashMap<String,String> ID2uniquePep = new HashMap<>();
        HashMap<String, HashSet<String>> ID2LRP = new HashMap<>();
        HashMap<String, HashSet<Integer>> ID2Positions = new HashMap<>();
        HashMap<String, HashSet<String>> ID2TopMods = new HashMap<>();
        HashMap<String, Double> ID2aScore = new HashMap<>();
        HashMap<String, String> ID2RawPep = new HashMap<>();
        HashMap<String, String> ID2Start = new HashMap<>();
        HashMap<String, String> ID2End = new HashMap<>();*/



        System.out.println("Complete");

        for (String ID: ids2Remove) {
            ID2Positions.remove(ID);
            ID2TopMods.remove(ID);
            ID2RawPep.remove(ID);
            ID2aScore.remove(ID);
            ID2Start.remove(ID);
            ID2End.remove(ID);
        }

        for (String lrp: ID2LRP.keySet()) {
            HashSet<String> ids = ID2LRP.get(lrp);
            for (String ID: ids2Remove) {
                if (ids.contains(ID)){
                    ids.remove(ID);
                    ID2LRP.put(lrp, ids);
                }
            }
        }




        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////  subsets & mapping

        System.out.println("Selecting UIDs, mapping to database, and generating neighbourhoods...");
        HashMap<Long, ArrayList<Integer>> numMapped = new HashMap<>();
        HashMap<Long,  ArrayList<Double>> avgSScore = new HashMap<>();
        HashMap<Long,  ArrayList<Double>> sumSScore = new HashMap<>();
        HashMap<Long,  Integer> nbhdIntegrated = new HashMap<>();
        HashMap<Long,  Integer> nbhdSize = new HashMap<>();
        HashMap<Long,  Integer> numUIDs = new HashMap<>();
        HashMap<Long,  ArrayList<Integer>> numUIDsMapped = new HashMap<>();
        HashMap<Long, String> nbhd2UID = new HashMap<>();
        HashMap<Long, String> nbhd2DispName = new HashMap<>();



        for (int j = 1; j < (repetitionNumber+1); j ++){
            System.out.println("\rCurrently performing replicate: " + j);

            //sample strategy
            // pick rand uids
            // make new ID2LRP
            // make new ID2TopMod
            // make new ID2Start
            // make new ID2End
            // make new ID2aScore

            Random random = new Random();
            HashSet<String> UIDs = new HashSet();

            while(UIDs.size() < subsetSize){
                int index = random.nextInt(ID2LRP.keySet().size());
                ArrayList<String> uidSet = new ArrayList<>(ID2LRP.keySet());
                UIDs.add(uidSet.get(index));
            }

            HashMap<String, HashSet<String>> ID2LRP_subset = new HashMap<>();
            HashMap<String, HashSet<String>> ID2TopMods_subset = new HashMap<>();
            HashMap<String, Double> ID2aScore_subset = new HashMap<>();
            HashMap<String, String> ID2Start_subset = new HashMap<>();
            HashMap<String, String> ID2End_subset = new HashMap<>();

            // pick rand num of peps from uid
            // no pick an index from a array with a distribution
            ArrayList<Long> distributionPercent = new ArrayList<>();
            ArrayList<Integer> distribution = new ArrayList<>();
            distributionPercent.add(Math.round(subsetSize*0.4));
            distributionPercent.add(Math.round(subsetSize*0.21));
            distributionPercent.add(Math.round(subsetSize*0.1));
            distributionPercent.add(Math.round(subsetSize*0.08));
            distributionPercent.add(Math.round(subsetSize*0.04));
            distributionPercent.add(Math.round(subsetSize*0.04));
            distributionPercent.add(Math.round(subsetSize*0.02));
            distributionPercent.add(Math.round(subsetSize*0.02));
            distributionPercent.add(Math.round(subsetSize*0.02));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));

            Integer count = 0;
            for (Long percent: distributionPercent) {
                count++;
                for (int i = 0; i < percent; i++) {
                    distribution.add(count);
                }
            }
            if(distribution.size() < subsetSize){
                Integer diff = subsetSize - distribution.size();
                for (int i = 0; i < diff ; i++) {
                    distribution.add(17);
                }
            }else if (distribution.size() > subsetSize){
                Integer diff = distribution.size() - subsetSize;
                for (int i = 0; i < diff ; i++) {
                    distribution.remove(distribution.size()-1);
                }
            }
            // now you have distribution =  [1,1,1,1,1,1,1,1,1,2,2,2,3,3,4,4,5,6] etc ..

            for(String randUID: UIDs){
                HashSet<String> pepIDs_subset = new HashSet<>();
                ArrayList<String> peptideIDs = new ArrayList<>(ID2LRP.get(randUID));
                int distIndex = random.nextInt(distribution.size());
                Integer numPeps = distribution.get(distIndex);

                if(peptideIDs.size() < numPeps){ // if less than the number of peptides available, take all
                    pepIDs_subset.addAll(peptideIDs);
                }else{ // if there are more peptides available than needed, pick random ones
                    // pick rand
                    int[] rand = new Random().ints(0, peptideIDs.size()).distinct().limit(numPeps).toArray();

                    for (int index: rand) {
                        pepIDs_subset.add(peptideIDs.get(index));

                    }
                }


                ID2LRP_subset.put(randUID, pepIDs_subset);

                for (String pepID: pepIDs_subset){
                    ID2TopMods_subset.put(pepID, ID2TopMods.get(pepID));
                    ID2aScore_subset.put(pepID, ID2aScore.get(pepID));
                    ID2Start_subset.put(pepID, ID2Start.get(pepID));
                    ID2End_subset.put(pepID, ID2End.get(pepID));
                }

            }

            // map to db
            mapToDB(graphDb,
                    ID2LRP_subset,
                    ID2TopMods_subset,
                    ID2Start_subset,
                    ID2End_subset,
                    ID2aScore_subset);



            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // generate neighbourhoods

            try (Transaction tx = graphDb.beginTx()) {


                // get the number of measured proteins and complexes
                int dbMeasured = 0;
                double numInDb = 0;
                ResourceIterator<Node> allProtandCmplxs = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
                for (ResourceIterator<Node> it = allProtandCmplxs; it.hasNext(); ) {
                    Node pe = it.next();
                    if(pe.hasProperty(PropertyType.SUPPORT_SCORE.toString())){
                        dbMeasured ++;
                    }
                    numInDb ++;
                }

                ResourceIterator<Node> allPeNodes = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
                for (ResourceIterator<Node> it = allPeNodes; it.hasNext(); ) {
                    Node physicalEntity = it.next();

                    // if pe is measured - gen nbhd
                    ArrayList<Double> nbhdSScores = new ArrayList<>();
                    Integer nbhdSizeInt = 0;
                    Integer nbhdMeasured = 0;
                    Integer numIntegrated = 0;
                    HashSet<String> nbhdUids = new HashSet<>();
                    HashSet<String> numUIDsMeasured = new HashSet<>();

                    if(physicalEntity.hasLabel(Label.label("Protein")) | physicalEntity.hasLabel(Label.label("Complex"))) {
                        HashSet<Node> nbhd = neighbourTraversal(physicalEntity, depth, graphDb);

                        // count the total neighbourhood size
                        // count number of measured things
                        for(Node node:nbhd){
                            if(node.hasLabel(Label.label("Protein")) | node.hasLabel(Label.label("Complex"))){
                                nbhdSizeInt ++;
                                //count integrated
                                if(node.hasProperty(PropertyType.INTEGRATED.toString())){
                                    numIntegrated++;
                                }
                            }

                            // get UID's from proteins and complexes in the network & count the # measured
                            if(node.hasLabel(Label.label("Protein"))){
                                Iterable<Relationship> uidRels = node.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                                for(Relationship uidRel: uidRels){
                                    nbhdUids.add(uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                                    if(uidRel.getStartNode().hasProperty(PropertyType.MAPPED.toString())){
                                        numUIDsMeasured.add(uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                                    }
                                }
                            }
                            if( node.hasLabel(Label.label("Complex"))){
                                String cplxUIDs = node.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                                String[] split = cplxUIDs.split(",\\s");
                                for (String uid: split) {
                                    Pattern p = Pattern.compile(UID_PATTERN);
                                    Matcher m = p.matcher(uid);
                                    if (m.find()) {
                                        nbhdUids.add(m.group(0));
                                    }
                                }
                                if(node.hasProperty(PropertyType.MAPPED.toString())){
                                    String mappedUIDs = node.getProperty(PropertyType.MAPPED.toString()).toString();
                                    String[] split1 = mappedUIDs.split(",\\s");
                                    for (String uid: split1){
                                        numUIDsMeasured.add(uid);
                                    }
                                }
                            }

                            // count measured entities
                            if (node.hasProperty(PropertyType.SUPPORT_SCORE.toString())){
                                if(Double.parseDouble(node.getProperty(PropertyType.SUPPORT_SCORE.toString()).toString()) > 0){
                                    nbhdMeasured ++;
                                }
                                nbhdSScores.add(Double.parseDouble(node.getProperty(PropertyType.SUPPORT_SCORE.toString()).toString()));
                            }
                        }

                        // add # measured per nbhd to hashmap
                        if(numMapped.containsKey(physicalEntity.getId())){
                            ArrayList<Integer> integers = numMapped.get(physicalEntity.getId());
                            integers.add(nbhdMeasured);
                            numMapped.put(physicalEntity.getId(), integers);
                        }else {
                            ArrayList<Integer> integers = new ArrayList<>();
                            integers.add(nbhdMeasured);
                            numMapped.put(physicalEntity.getId(), integers);
                        }

                        // add # measured per nbhd to hashmap
                        if(numUIDsMapped.containsKey(physicalEntity.getId())){
                            ArrayList<Integer> integers = numUIDsMapped.get(physicalEntity.getId());
                            integers.add(numUIDsMeasured.size());
                            numUIDsMapped.put(physicalEntity.getId(), integers);
                        }else {
                            ArrayList<Integer> integers = new ArrayList<>();
                            integers.add(numUIDsMeasured.size());
                            numUIDsMapped.put(physicalEntity.getId(), integers);
                        }

                        nbhdSize.put(physicalEntity.getId(), nbhdSizeInt);

                        nbhdIntegrated.put(physicalEntity.getId(), numIntegrated);

                        numUIDs.put(physicalEntity.getId(), nbhdUids.size());


                        // AVG & SUM & STDDEV
                        Double ssAvg = 0.0;
                        Double ssSum = 0.0;
                        for(Integer i = 0; i < nbhdSScores.size(); i++){
                            if(nbhdSScores.get(i) >= 0){
                                ssAvg += nbhdSScores.get(i);
                                ssSum += nbhdSScores.get(i);
                            }
                        }
                        if(!ssSum.equals(0.0)){
                            ssAvg = ssAvg/nbhdSScores.size();
                        }

                        // add avg SScore per nbhd to hashmap
                        ssAvg = (Math.round(ssAvg*100.00))/100.00;
                        if(avgSScore.containsKey(physicalEntity.getId())){
                            ArrayList<Double> doubles = avgSScore.get(physicalEntity.getId());
                            doubles.add(ssAvg);
                            avgSScore.put(physicalEntity.getId(), doubles);
                        }else {
                            ArrayList<Double> doubles = new ArrayList<>();
                            doubles.add(ssAvg);
                            avgSScore.put(physicalEntity.getId(), doubles);
                        }

                        ssSum = (Math.round(ssSum*100.00))/100.00;
                        if(sumSScore.containsKey(physicalEntity.getId())){
                            ArrayList<Double> doubles = sumSScore.get(physicalEntity.getId());
                            doubles.add(ssSum);
                            sumSScore.put(physicalEntity.getId(), doubles);
                        }else {
                            ArrayList<Double> doubles = new ArrayList<>();
                            doubles.add(ssSum);
                            sumSScore.put(physicalEntity.getId(), doubles);
                        }

                        ////////////// gathering outputs
                        String uid = "";

                        if(physicalEntity.hasProperty(PropertyType.SCORED_BY.toString())){
                            uid = physicalEntity.getProperty(PropertyType.SCORED_BY.toString()).toString();
                        } else if (physicalEntity.hasProperty(PropertyType.UNIPROT_ID.toString())){
                            uid = physicalEntity.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                        }else{
                            Iterable<Relationship> uidRels = physicalEntity.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                            if(getLength(uidRels) == 1){
                                for(Relationship uidRel: uidRels){
                                    uid = uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                                }
                            }else{
                                for(Relationship uidRel: uidRels){
                                    uid = uid + uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString() + ", ";
                                }
                            }
                        }
                        nbhd2UID.put(physicalEntity.getId(), uid);
                        nbhd2DispName.put(physicalEntity.getId(), physicalEntity.getProperty(PropertyType.DISPLAY_NAME.toString()).toString());
                    }
                }

                // reset scores
                ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
                for (ResourceIterator<Node> it = proteins; it.hasNext(); ) {
                    Node protein = it.next();
                    protein.removeProperty(PropertyType.SUPPORT_SCORE.toString());
                    protein.removeProperty(PropertyType.ABUNDANCE_SCORE.toString());
                    protein.removeProperty(PropertyType.SCORED_BY.toString());
                }

                ResourceIterator<Node> complexs = graphDb.findNodes(Label.label("Complex"));
                for (ResourceIterator<Node> it = complexs; it.hasNext(); ) {
                    Node complex = it.next();
                    complex.removeProperty(PropertyType.SUPPORT_SCORE.toString());
                    complex.removeProperty(PropertyType.ABUNDANCE_SCORE.toString());
                    complex.removeProperty(PropertyType.SCORED_BY.toString());
                    complex.removeProperty(PropertyType.MAPPED.toString());
                }

                ResourceIterator<Node> uids = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
                for (ResourceIterator<Node> it = uids; it.hasNext(); ) {
                    Node uid = it.next();
                    uid.removeProperty(PropertyType.MAPPED.toString());
                }

                tx.success();
            }
        }

        // get Uniprot accesison, position,  Raw Peptide, log2 ratio
        // clear file if already exits
        FileWriter fstream0 = new FileWriter(outputFile + "/EmpiricalDist.tsv");
        BufferedWriter out0 = new BufferedWriter(fstream0);
        out0.write("nbhdBaseDBID\tnbhdBaseUID\tnbhdBaseDisplayName\tnumEntities\tNumIntegrated\tnumUIDs\tnumUIDsMapped\tnbhdMeasured\tavgSScore\tsumSScore\n" );
        out0.close();
        //Create File to write inputs/Outputs to
        FileWriter fstream = new FileWriter(outputFile + "/EmpiricalDist.tsv", true);
        BufferedWriter out = new BufferedWriter(fstream);

        for (Long nbhd: numMapped.keySet()) {
            out.write(nbhd +"\t"
                    + nbhd2UID.get(nbhd)
                    +"\t" + nbhd2DispName.get(nbhd)
                    +"\t" + nbhdSize.get(nbhd)
                    +"\t" + nbhdIntegrated.get(nbhd)
                    +"\t" + numUIDs.get(nbhd)
                    +"\t" + numUIDsMapped.get(nbhd)
                    +"\t" + numMapped.get(nbhd)
                    +"\t" + avgSScore.get(nbhd)
                    +"\t" + sumSScore.get(nbhd)
                    +"\n");

        }

        out.close();
        graphDb.shutdown();
    }

    void mouseEmpiricalNullDistribution(File mouseDBFile, Integer depth, Integer subsetSize, Integer repetitionNumber ) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        // reset scores
        try (Transaction tx = graphDb.beginTx()) {
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            for (ResourceIterator<Node> it = proteins; it.hasNext(); ) {
                Node protein = it.next();
                protein.removeProperty(PropertyType.SUPPORT_SCORE.toString());
                protein.removeProperty(PropertyType.ABUNDANCE_SCORE.toString());
                protein.removeProperty(PropertyType.SCORED_BY.toString());
            }

            ResourceIterator<Node> complexs = graphDb.findNodes(Label.label("Complex"));
            for (ResourceIterator<Node> it = complexs; it.hasNext(); ) {
                Node complex = it.next();
                complex.removeProperty(PropertyType.SUPPORT_SCORE.toString());
                complex.removeProperty(PropertyType.ABUNDANCE_SCORE.toString());
                complex.removeProperty(PropertyType.SCORED_BY.toString());
                complex.removeProperty(PropertyType.MAPPED.toString());
            }

            ResourceIterator<Node> uids = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            for (ResourceIterator<Node> it = uids; it.hasNext(); ) {
                Node uid = it.next();
                uid.removeProperty(PropertyType.MAPPED.toString());
            }
            tx.success();
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Reading file in and getting hashmaps for mapping

        // Read in data into Dict
        BufferedReader BR = new BufferedReader(new FileReader(mouseDBFile));
        String line = "";

        HashMap<String, HashSet<String>> ID2LRP = new HashMap<>();
        HashMap<String, Double> ID2aScore = new HashMap<>();
        HashMap<String, HashSet<String>> ID2TopMods = new HashMap<>();
        HashMap<String, String> ID2Start = new HashMap<>();
        HashMap<String, String> ID2End = new HashMap<>();
        HashSet<String> ids2Remove = new HashSet<>();


        Integer leadingProtCol =0;
        Integer aScoreCol = 0;
        Integer modPepCol = 0;

        System.out.println("Reading in all qPhos Data");

        String species = getSpecies(graphDb);
        HashMap<String, String> fastaMap = new HashMap<>();
        if(species.equalsIgnoreCase("human")){
            // read in human fasta
            fastaMap = fasta2dict(humanUniProt);
        }else if(species.equalsIgnoreCase("mouse")){
            // read in mouse fasta
            fastaMap = fasta2dict(mouseUniProt);
        }

        Integer countID = 0;
        while ((line = BR.readLine()) != null) {
            String[] columns = line.split("\t");
            Integer colLen = columns.length;
            if(countID == 0){ // if first line in file
                for (int i = 0; i < colLen; i++) {
                    if (columns[i].equalsIgnoreCase("LRP")){
                        leadingProtCol = i;
                    }
                    if(columns[i].equals("Log2INT")){
                        aScoreCol = i;
                    }
                    if(columns[i].equals("mod.pep.seq")){
                        modPepCol = i;
                    }
                }

                countID++;
            }else{ //for all other lines
                // extract LRP's
                if(ID2LRP.containsKey(columns[leadingProtCol])){
                    HashSet<String> temp = ID2LRP.get(columns[leadingProtCol]);
                    temp.add(String.valueOf(countID));
                    ID2LRP.put(columns[leadingProtCol], temp);
                }else{
                    HashSet<String> temp = new HashSet<>();
                    temp.add(String.valueOf(countID));
                    ID2LRP.put(columns[leadingProtCol], temp);
                }

                // extract aScore (log2Ints)
                ID2aScore.put(String.valueOf(countID), Double.valueOf(columns[aScoreCol]));

                //generate Start
                String start = generateStartOfPep(fastaMap, columns[leadingProtCol],columns[modPepCol]);
                if(start  == ""){
                    ids2Remove.add(countID.toString());
                }else{
                    ID2Start.put(String.valueOf(countID), start);

                    // generate End
                    String end = generateEndOfPep(start, columns[modPepCol]);
                    ID2End.put(String.valueOf(countID), end);

                    // generate Top Mods
                    HashSet<String> mods = generateModList(countID,columns[modPepCol], ID2Start);
                    ID2TopMods.put(String.valueOf(countID), mods);
                }

                countID++;
            }
        }

        System.out.println("Complete");

        for (String ID: ids2Remove) {
            ID2TopMods.remove(ID);
            ID2aScore.remove(ID);
            ID2Start.remove(ID);
            ID2End.remove(ID);
        }

        for (String lrp: ID2LRP.keySet()) {
            HashSet<String> ids = ID2LRP.get(lrp);
            for (String ID: ids2Remove) {
                if (ids.contains(ID)){
                    ids.remove(ID);
                    ID2LRP.put(lrp, ids);
                }
            }
        }




        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////  subsets & mapping

        System.out.println("Selecting UIDs, mapping to database, and generating neighbourhoods...");
        HashMap<Long, ArrayList<Integer>> numMapped = new HashMap<>();
        HashMap<Long,  ArrayList<Double>> avgSScore = new HashMap<>();
        HashMap<Long,  ArrayList<Double>> sumSScore = new HashMap<>();
        HashMap<Long,  Integer> nbhdIntegrated = new HashMap<>();
        HashMap<Long,  Integer> nbhdSize = new HashMap<>();
        HashMap<Long,  Integer> numUIDs = new HashMap<>();
        HashMap<Long,  ArrayList<Integer>> numUIDsMapped = new HashMap<>();
        HashMap<Long, String> nbhd2UID = new HashMap<>();
        HashMap<Long, String> nbhd2DispName = new HashMap<>();

        for (int j = 1; j < (repetitionNumber+1); j ++){
            System.out.println("\rCurrently performing replicate: " + j);

            //sample strategy
            // pick rand uids
            // make new ID2LRP
            // make new ID2TopMod
            // make new ID2Start
            // make new ID2End
            // make new ID2aScore

            Random random = new Random();
            HashSet<String> UIDs = new HashSet();

            while(UIDs.size() < subsetSize){
                int index = random.nextInt(ID2LRP.keySet().size());
                ArrayList<String> uidSet = new ArrayList<>(ID2LRP.keySet());
                UIDs.add(uidSet.get(index));
            }

            HashMap<String, HashSet<String>> ID2LRP_subset = new HashMap<>();
            HashMap<String, HashSet<String>> ID2TopMods_subset = new HashMap<>();
            HashMap<String, Double> ID2aScore_subset = new HashMap<>();
            HashMap<String, String> ID2Start_subset = new HashMap<>();
            HashMap<String, String> ID2End_subset = new HashMap<>();

            // pick rand num of peps from uid
            // no pick an index from a array with a distribution
            ArrayList<Long> distributionPercent = new ArrayList<>();
            ArrayList<Integer> distribution = new ArrayList<>();
            distributionPercent.add(Math.round(subsetSize*0.4));
            distributionPercent.add(Math.round(subsetSize*0.21));
            distributionPercent.add(Math.round(subsetSize*0.1));
            distributionPercent.add(Math.round(subsetSize*0.08));
            distributionPercent.add(Math.round(subsetSize*0.04));
            distributionPercent.add(Math.round(subsetSize*0.04));
            distributionPercent.add(Math.round(subsetSize*0.02));
            distributionPercent.add(Math.round(subsetSize*0.02));
            distributionPercent.add(Math.round(subsetSize*0.02));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));
            distributionPercent.add(Math.round(subsetSize*0.01));

            Integer count = 0;
            for (Long percent: distributionPercent) {
                count++;
                for (int i = 0; i < percent; i++) {
                    distribution.add(count);
                }
            }
            if(distribution.size() < subsetSize){
                Integer diff = subsetSize - distribution.size();
                for (int i = 0; i < diff ; i++) {
                    distribution.add(17);
                }
            }else if (distribution.size() > subsetSize){
                Integer diff = distribution.size() - subsetSize;
                for (int i = 0; i < diff ; i++) {
                    distribution.remove(distribution.size()-1);
                }
            }
            // now you have distribution =  [1,1,1,1,1,1,1,1,1,2,2,2,3,3,4,4,5,6] etc ..

            for(String randUID: UIDs){
                HashSet<String> pepIDs_subset = new HashSet<>();
                ArrayList<String> peptideIDs = new ArrayList<>(ID2LRP.get(randUID));
                int distIndex = random.nextInt(distribution.size());
                Integer numPeps = distribution.get(distIndex);

                if(peptideIDs.size() < numPeps){ // if less than the number of peptides available, take all
                    pepIDs_subset.addAll(peptideIDs);
                }else{ // if there are more peptides available than needed, pick random ones
                    // pick rand
                    int[] rand = new Random().ints(0, peptideIDs.size()).distinct().limit(numPeps).toArray();

                    for (int index: rand) {
                        pepIDs_subset.add(peptideIDs.get(index));

                    }
                }


                ID2LRP_subset.put(randUID, pepIDs_subset);

                for (String pepID: pepIDs_subset){
                    ID2TopMods_subset.put(pepID, ID2TopMods.get(pepID));
                    ID2aScore_subset.put(pepID, ID2aScore.get(pepID));
                    ID2Start_subset.put(pepID, ID2Start.get(pepID));
                    ID2End_subset.put(pepID, ID2End.get(pepID));
                }

            }


            // map to db
            mapToDB(graphDb,
                    ID2LRP_subset,
                    ID2TopMods_subset,
                    ID2Start_subset,
                    ID2End_subset,
                    ID2aScore_subset);



            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // generate neighbourhoods

            try (Transaction tx = graphDb.beginTx()) {


                // get the number of measured proteins and complexes
                int dbMeasured = 0;
                double numInDb = 0;
                ResourceIterator<Node> allProtandCmplxs = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
                for (ResourceIterator<Node> it = allProtandCmplxs; it.hasNext(); ) {
                    Node pe = it.next();
                    if(pe.hasProperty(PropertyType.SUPPORT_SCORE.toString())){
                        dbMeasured ++;
                    }
                    numInDb ++;
                }

                ResourceIterator<Node> allPeNodes = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
                for (ResourceIterator<Node> it = allPeNodes; it.hasNext(); ) {
                    Node physicalEntity = it.next();

                    // if pe is measured - gen nbhd
                    ArrayList<Double> nbhdSScores = new ArrayList<>();
                    Integer nbhdSizeInt = 0;
                    Integer nbhdMeasured = 0;
                    Integer numIntegrated = 0;
                    HashSet<String> nbhdUids = new HashSet<>();
                    HashSet<String> numUIDsMeasured = new HashSet<>();

                    if(physicalEntity.hasLabel(Label.label("Protein")) | physicalEntity.hasLabel(Label.label("Complex"))) {
                        HashSet<Node> nbhd = neighbourTraversal(physicalEntity, depth, graphDb);

                        // count the total neighbourhood size
                        // count number of measured things
                        for(Node node:nbhd){
                            if(node.hasLabel(Label.label("Protein")) | node.hasLabel(Label.label("Complex"))){
                                nbhdSizeInt ++;
                                //count integrated
                                if(node.hasProperty(PropertyType.INTEGRATED.toString())){
                                    numIntegrated++;
                                }
                            }

                            // get UID's from proteins and complexes in the network & count the # measured
                            if(node.hasLabel(Label.label("Protein"))){
                                Iterable<Relationship> uidRels = node.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                                for(Relationship uidRel: uidRels){
                                    nbhdUids.add(uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                                    if(uidRel.getStartNode().hasProperty(PropertyType.MAPPED.toString())){
                                        numUIDsMeasured.add(uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                                    }
                                }
                            }
                            if( node.hasLabel(Label.label("Complex"))){
                                String cplxUIDs = node.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                                String[] split = cplxUIDs.split(",\\s");
                                for (String uid: split) {
                                    Pattern p = Pattern.compile(UID_PATTERN);
                                    Matcher m = p.matcher(uid);
                                    if (m.find()) {
                                        nbhdUids.add(m.group(0));
                                    }
                                }
                                if(node.hasProperty(PropertyType.MAPPED.toString())){
                                    String mappedUIDs = node.getProperty(PropertyType.MAPPED.toString()).toString();
                                    String[] split1 = mappedUIDs.split(",\\s");
                                    for (String uid: split1){
                                        numUIDsMeasured.add(uid);
                                    }
                                }
                            }

                            // count measured entities
                            if (node.hasProperty(PropertyType.SUPPORT_SCORE.toString())){
                                if(Double.parseDouble(node.getProperty(PropertyType.SUPPORT_SCORE.toString()).toString()) > 0){
                                    nbhdMeasured ++;
                                }
                                nbhdSScores.add(Double.parseDouble(node.getProperty(PropertyType.SUPPORT_SCORE.toString()).toString()));
                            }
                        }

                        // add # measured per nbhd to hashmap
                        if(numMapped.containsKey(physicalEntity.getId())){
                            ArrayList<Integer> integers = numMapped.get(physicalEntity.getId());
                            integers.add(nbhdMeasured);
                            numMapped.put(physicalEntity.getId(), integers);
                        }else {
                            ArrayList<Integer> integers = new ArrayList<>();
                            integers.add(nbhdMeasured);
                            numMapped.put(physicalEntity.getId(), integers);
                        }

                        // add # measured per nbhd to hashmap
                        if(numUIDsMapped.containsKey(physicalEntity.getId())){
                            ArrayList<Integer> integers = numUIDsMapped.get(physicalEntity.getId());
                            integers.add(numUIDsMeasured.size());
                            numUIDsMapped.put(physicalEntity.getId(), integers);
                        }else {
                            ArrayList<Integer> integers = new ArrayList<>();
                            integers.add(numUIDsMeasured.size());
                            numUIDsMapped.put(physicalEntity.getId(), integers);
                        }

                        nbhdSize.put(physicalEntity.getId(), nbhdSizeInt);

                        nbhdIntegrated.put(physicalEntity.getId(), numIntegrated);

                        numUIDs.put(physicalEntity.getId(), nbhdUids.size());


                        // AVG & SUM & STDDEV
                        Double ssAvg = 0.0;
                        Double ssSum = 0.0;
                        for(Integer i = 0; i < nbhdSScores.size(); i++){
                            if(nbhdSScores.get(i) >= 0){
                                ssAvg += nbhdSScores.get(i);
                                ssSum += nbhdSScores.get(i);
                            }
                        }
                        if(!ssSum.equals(0.0)){
                            ssAvg = ssAvg/nbhdSScores.size();
                        }

                        // add avg SScore per nbhd to hashmap
                        ssAvg = (Math.round(ssAvg*100.00))/100.00;
                        if(avgSScore.containsKey(physicalEntity.getId())){
                            ArrayList<Double> doubles = avgSScore.get(physicalEntity.getId());
                            doubles.add(ssAvg);
                            avgSScore.put(physicalEntity.getId(), doubles);
                        }else {
                            ArrayList<Double> doubles = new ArrayList<>();
                            doubles.add(ssAvg);
                            avgSScore.put(physicalEntity.getId(), doubles);
                        }

                        ssSum = (Math.round(ssSum*100.00))/100.00;
                        if(sumSScore.containsKey(physicalEntity.getId())){
                            ArrayList<Double> doubles = sumSScore.get(physicalEntity.getId());
                            doubles.add(ssSum);
                            sumSScore.put(physicalEntity.getId(), doubles);
                        }else {
                            ArrayList<Double> doubles = new ArrayList<>();
                            doubles.add(ssSum);
                            sumSScore.put(physicalEntity.getId(), doubles);
                        }

                        ////////////// gathering outputs
                        String uid = "";

                        if(physicalEntity.hasProperty(PropertyType.SCORED_BY.toString())){
                            uid = physicalEntity.getProperty(PropertyType.SCORED_BY.toString()).toString();
                        } else if (physicalEntity.hasProperty(PropertyType.UNIPROT_ID.toString())){
                            uid = physicalEntity.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                        }else{
                            Iterable<Relationship> uidRels = physicalEntity.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                            if(getLength(uidRels) == 1){
                                for(Relationship uidRel: uidRels){
                                    uid = uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                                }
                            }else{
                                for(Relationship uidRel: uidRels){
                                    uid = uid + uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString() + ", ";
                                }
                            }
                        }
                        nbhd2UID.put(physicalEntity.getId(), uid);
                        nbhd2DispName.put(physicalEntity.getId(), physicalEntity.getProperty(PropertyType.DISPLAY_NAME.toString()).toString());
                    }
                }

                // reset scores
                ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
                for (ResourceIterator<Node> it = proteins; it.hasNext(); ) {
                    Node protein = it.next();
                    protein.removeProperty(PropertyType.SUPPORT_SCORE.toString());
                    protein.removeProperty(PropertyType.ABUNDANCE_SCORE.toString());
                    protein.removeProperty(PropertyType.SCORED_BY.toString());
                }

                ResourceIterator<Node> complexs = graphDb.findNodes(Label.label("Complex"));
                for (ResourceIterator<Node> it = complexs; it.hasNext(); ) {
                    Node complex = it.next();
                    complex.removeProperty(PropertyType.SUPPORT_SCORE.toString());
                    complex.removeProperty(PropertyType.ABUNDANCE_SCORE.toString());
                    complex.removeProperty(PropertyType.SCORED_BY.toString());
                    complex.removeProperty(PropertyType.MAPPED.toString());
                }

                ResourceIterator<Node> uids = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
                for (ResourceIterator<Node> it = uids; it.hasNext(); ) {
                    Node uid = it.next();
                    uid.removeProperty(PropertyType.MAPPED.toString());
                }

                tx.success();
            }
        }


        // get Uniprot accesison, position,  Raw Peptide, log2 ratio
        // clear file if already exits
        FileWriter fstream0 = new FileWriter(outputFile + "/EmpiricalDist.tsv");
        BufferedWriter out0 = new BufferedWriter(fstream0);
        out0.write("nbhdBaseDBID\tnbhdBaseUID\tnbhdBaseDisplayName\tnumEntities\tNumIntegrated\tnumUIDs\tnumUIDsMapped\tnbhdMeasured\tavgSScore\tsumSScore\n" );
        out0.close();
        //Create File to write inputs/Outputs to
        FileWriter fstream = new FileWriter(outputFile + "/EmpiricalDist.tsv", true);
        BufferedWriter out = new BufferedWriter(fstream);

        for (Long nbhd: numMapped.keySet()) {
            out.write(nbhd +"\t"
                    + nbhd2UID.get(nbhd)
                    +"\t" + nbhd2DispName.get(nbhd)
                    +"\t" + nbhdSize.get(nbhd)
                    +"\t" + nbhdIntegrated.get(nbhd)
                    +"\t" + numUIDs.get(nbhd)
                    +"\t" + numUIDsMapped.get(nbhd)
                    +"\t" + numMapped.get(nbhd)
                    +"\t" + avgSScore.get(nbhd)
                    +"\t" + sumSScore.get(nbhd)
                    +"\n");

        }

        out.close();
        graphDb.shutdown();

    }

    private void mapToDB(GraphDatabaseService graphDb,
                         HashMap<String, HashSet<String>> ID2LRP,
                         HashMap<String, HashSet<String>> ID2TopMods,
                         HashMap<String, String> ID2Start,
                         HashMap<String, String> ID2End,
                         HashMap<String, Double> ID2aScore){

        try (Transaction tx = graphDb.beginTx()) {

            Integer proteoformsScored = 0;
            double pepsSpanningModTemp = 0;
            double pepsSpanningMod = 0;
            Integer pepsMatched = 0;
            Integer UIDsMatched = 0;
            Integer UIDsNotMatched = 0;
            Integer PepsNotMatched = 0;

            // for each LRP in the data

            for(String LRP: ID2LRP.keySet()){
                HashSet<String> PeptideIDs = ID2LRP.get(LRP);

                // get the UID from the DB
                Node UID = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), LRP);


                if (UID == null) {
                    UIDsNotMatched ++;
                    PepsNotMatched += PeptideIDs.size();
                    continue;
                } else {
                    UIDsMatched ++;
                    pepsMatched += PeptideIDs.size();

                    UID.setProperty(PropertyType.MAPPED.toString(), "Mapped");

                    // GET UID
                    Iterable<Relationship> protRelationships = UID.getRelationships(RelTypes.ID_BELONGS_TO);

                    Integer count = 0;
                    HashSet<String> UIDmodList = new HashSet<>();
                    HashMap<Integer, String> uidModList = new HashMap<>();

                    // GET POIs
                    for (Relationship protRltnshp : protRelationships) {
                        Node protein = protRltnshp.getEndNode();
                        Iterable<Relationship> phosRelationships = protein.getRelationships(RelTypes.PHOSPHORYLATION);
                        for (Relationship phosRltnshp : phosRelationships) {
                            // GET ALL PHOSPHOTYLATIONS
                            Node phos = phosRltnshp.getStartNode();
                            String templocation = phos.getProperty(PropertyType.LOCATION.toString()).toString();
                            if (!templocation.equals("unknown") & !templocation.equals("NA")){
                                Integer location = Integer.valueOf(phos.getProperty(PropertyType.LOCATION.toString()).toString());
                                String type = phos.getProperty(PropertyType.TYPE.toString()).toString().replaceAll("p\\_", "");
                                String rxmMod = type + "_" + location;
                                UIDmodList.add(rxmMod);
                                uidModList.put(location, type);
                            }
                        }
                    }


                    // GET PF MODs
                    for (Relationship protRltnshp : protRelationships) {
                        double proteoformScore = 0.0;
                        Node protein = protRltnshp.getEndNode();
                        Iterable<Relationship> phosRelationships = protein.getRelationships(RelTypes.PHOSPHORYLATION);
                        HashSet<String> proteoformModList = new HashSet<>();
                        HashMap<Integer, String> pfModList = new HashMap<>();
                        for (Relationship phosRltnshp : phosRelationships) {
                            // GET ALL PHOSPHOTYLATIONS
                            Node phos = phosRltnshp.getStartNode();
                            String templocation = phos.getProperty(PropertyType.LOCATION.toString()).toString();

                            if (!templocation.equals("unknown") & !templocation.equals("NA")){
                                Integer location = Integer.valueOf(phos.getProperty(PropertyType.LOCATION.toString()).toString());
                                String type = phos.getProperty(PropertyType.TYPE.toString()).toString().replaceAll("p\\_", "");
                                String rxmMod = type + "_" + location;
                                proteoformModList.add(rxmMod);
                                pfModList.put(location, type);
                            }

                        }



                        HashMap<String, Double> pepID2sScore = new HashMap<>();
                        for(String PeptideID: PeptideIDs) {
                            //System.out.println(PeptideID + "\t" + ID2LRP.get(LRP) + "\t" + ID2TopMods.get(PeptideID) + "\t" +ID2Start.get(PeptideID));
                            HashSet<String> pepMods = ID2TopMods.get(PeptideID);
                            Integer start = Integer.valueOf(ID2Start.get(PeptideID));
                            Integer end = Integer.valueOf(ID2End.get(PeptideID));


                            // get mods covered by peptide
                            HashSet<String> POIwithinPep = new HashSet<>();
                            for(Integer uidLocation: uidModList.keySet()){
                                if(uidLocation >= start & uidLocation <= end){
                                    POIwithinPep.add(uidModList.get(uidLocation) + "_"+ uidLocation);
                                }
                            }

                            if(!POIwithinPep.isEmpty()){
                                pepsSpanningModTemp ++;
                            }


                            double unModPeptideScore = 0.0;
                            double peptideScore = 0.0;
                            if(POIwithinPep.isEmpty()){ // no mods on uid
                                if(!pepMods.isEmpty()){ // if mods in pep
                                    peptideScore+= 0.9;
                                }else{ // no mods in pep (perf score without multiplying by 1.5)
                                    unModPeptideScore += 1.0;
                                }
                            }else{

                                for(String POI: POIwithinPep){ //
                                    if((proteoformModList.contains(POI) & pepMods.contains(POI)) // match
                                            | (!proteoformModList.contains(POI) & !pepMods.contains(POI))){
                                        peptideScore += 1.0;
                                        //System.out.println("Modified pf and mods in pep match += 1");
                                    }else if((proteoformModList.contains(POI) & !pepMods.contains(POI)) // mismatch
                                            | (!proteoformModList.contains(POI) & pepMods.contains(POI))){
                                        peptideScore += 0.5;
                                        //System.out.println("Modified pf but pep mods dont match += 0.5");
                                    }
                                    //System.out.print(peptideScore +" + ");
                                }


                                peptideScore = peptideScore/ POIwithinPep.size();
                                if(peptideScore == 1){ // if an exact match multiply by 1.5
                                    if(!pepMods.isEmpty()){
                                        peptideScore = peptideScore * 1.5;
                                        //System.out.println("Modified pep had a perfect score so *1.5");
                                    }
                                }

                                // make a copy of all mods in a peptide
                                HashSet<String> pepModCopy = new HashSet<>();
                                for(String pepMod: pepMods){
                                    pepModCopy.add(pepMod);
                                }

                                pepModCopy.removeAll(UIDmodList);

                                if(!pepModCopy.isEmpty()){ // if pep has mods that are not POI, not exact match changed * 0.9
                                    peptideScore = peptideScore*0.9;
                                }
                            }
                            //System.out.println(" = " + peptideScore);
                            //System.out.println("um = " + unModPeptideScore);

                            pepID2sScore.put(PeptideID, (peptideScore+unModPeptideScore));
                            proteoformScore += peptideScore;
                            proteoformScore += unModPeptideScore;
                        }
                        proteoformsScored ++;
                        proteoformScore = proteoformScore/PeptideIDs.size();

                        //System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()) +" , "+pepID2sScore);

                        // If protein already has a score - must have multi UIDs
                        // so take the uid that gives it the highest score
                        //System.out.println( "PF score: "+ proteoformScore);
                        if(protein.hasProperty(PropertyType.SUPPORT_SCORE.toString())){
                            double currentScore = Double.parseDouble(protein.getProperty(PropertyType.SUPPORT_SCORE.toString()).toString());
                            if(proteoformScore > currentScore){
                                protein.setProperty(PropertyType.SUPPORT_SCORE.toString(), proteoformScore);
                                protein.setProperty(PropertyType.SCORED_BY.toString(), LRP);
                            }
                        }else{
                            protein.setProperty(PropertyType.SUPPORT_SCORE.toString(), proteoformScore);
                        }


                        pepsSpanningModTemp = pepsSpanningModTemp/getLength(protRelationships);
                        pepsSpanningMod += pepsSpanningModTemp;
                        pepsSpanningModTemp =0;

                        // Setting abundance score: // TODO ??
                        // find the Peptide ID's with the biggest score
                        Map.Entry<String, Double> maxEntry = null;
                        for (Map.Entry<String, Double> entry: pepID2sScore.entrySet()){
                            if(maxEntry == null || entry.getValue().compareTo(maxEntry.getValue()) > 0 ){
                                maxEntry = entry;
                            }
                        }

                        // may be multi pep's with largest score, add all to hashset
                        HashSet<String> maxPeps = new HashSet<>();
                        for(String pep: pepID2sScore.keySet()){
                            Double score = pepID2sScore.get(pep);
                            if (score.equals(maxEntry.getValue())){
                                maxPeps.add(pep);
                            }
                        }

                        // if multi get avg e_score of highest ranking peps (per pf)
                        // else take highest e_score
                        Double eScore = 0.0;
                        if(maxPeps.size() >1){ // if multi high supported - get avg

                            // create array of eScores to avg
                            double eScores[] = new double[maxPeps.size()];
                            int pepcount = 0;
                            for(String pep: maxPeps){
                                eScores[pepcount] = ID2aScore.get(pep);
                                pepcount ++;
                            }
                            // avg escores
                            double sum = 0;
                            for(double i: eScores){
                                sum += i;
                            }
                            eScore = sum/eScores.length;

                        }else{ // maxPeps.size == 1
                            for(String pep: maxPeps){ // should only be one
                                eScore = ID2aScore.get(pep);
                            }
                        }
                        protein.setProperty(PropertyType.ABUNDANCE_SCORE.toString(), eScore);
                    }
                }
            }

            tx.success();
            mapComplexs(graphDb, "SUPPORT_SCORE", "ABUNDANCE_SCORE",  "SCORED_BY", "MAPPED");

        }

    }

    void nbhdAnalysis(Integer depth, String experiment, File empiricalDistribution) throws IOException {

        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        ////////////////////////////////////////////// Reading in empirical dists //////////////////////////////////////////////
        BufferedReader BR = new BufferedReader(new FileReader(empiricalDistribution));
        String line = "";

        Integer nbhdBaseDBIDCol =0;
        Integer numEntitiesCol = 0;
        Integer numIntegratedCol = 0;
        Integer numUIDsCol = 0;
        Integer numUIDsMappedCol = 0;
        Integer nbhdMesuredCol = 0;
        Integer avgSScoreCol = 0;
        Integer sumSScoreCol = 0;

        System.out.println("Reading in all distributions per Neighbourhood");

        Integer countID = 0;

        HashMap<Long, ArrayList<Integer>> numUIDSMappedQPhos = new HashMap<>();
        HashMap<Long, ArrayList<Integer>> numEntitiesMappedQPhos = new HashMap<>();
        HashMap<Long, ArrayList<Double>> avgSScoreQPhos = new HashMap<>();
        HashMap<Long, ArrayList<Double>> sumSScoreQPhos = new HashMap<>();

        while ((line = BR.readLine()) != null) {
            String[] columns = line.split("\t");
            Integer colLen = columns.length;
            if (countID == 0) { // if first line in file
                for (int i = 0; i < colLen; i++) {
                    if (columns[i].equalsIgnoreCase("nbhdBaseDBID")) {
                        nbhdBaseDBIDCol = i;
                    }
                    if (columns[i].equals("numEntities")) {
                        numEntitiesCol = i;
                    }
                    if (columns[i].equals("NumIntegrated")) {
                        numIntegratedCol = i;
                    }
                    if (columns[i].equals("numUIDs")) {
                        numUIDsCol = i;
                    }
                    if (columns[i].equals("numUIDsMapped")) {
                        numUIDsMappedCol = i;
                    }
                    if (columns[i].equals("nbhdMeasured")) {
                        nbhdMesuredCol = i;
                    }
                    if (columns[i].equals("avgSScore")) {
                        avgSScoreCol = i;
                    }
                    if (columns[i].equals("sumSScore")) {
                        sumSScoreCol = i;
                    }

                }
                countID++;
            } else { //for all other lines

                String nbhdBaseDBID = columns[nbhdBaseDBIDCol];
                String numEntities = columns[numEntitiesCol];
                String numIntegrated = columns[numIntegratedCol];
                String numUIDs = columns[numUIDsCol];
                String numUIDsMapped = columns[numUIDsMappedCol];
                String nbhdMesured = columns[nbhdMesuredCol];
                String avgSScore = columns[avgSScoreCol];
                String sumSScore = columns[sumSScoreCol];

                numUIDsMapped = numUIDsMapped.replace("[", "");
                numUIDsMapped = numUIDsMapped.replace("]", "");
                ArrayList<Integer>  numUIDsMappedSet = new ArrayList<>();
                for (String s :  numUIDsMapped.split(",\\s") ){
                    numUIDsMappedSet.add(Integer.valueOf(s));
                }

                ArrayList<Integer>  nbhdMeasuredSet = new ArrayList<>();
                nbhdMesured = nbhdMesured.replace("[", "");
                nbhdMesured = nbhdMesured.replace("]", "");
                for (String s :  nbhdMesured.split(",\\s") ){
                    nbhdMeasuredSet.add(Integer.valueOf(s));
                }

                ArrayList<Double>  avgSScoreSet = new ArrayList<>();
                avgSScore = avgSScore.replace("[", "");
                avgSScore = avgSScore.replace("]", "");
                for (String s :  avgSScore.split(",\\s") ){
                    avgSScoreSet.add(Double.valueOf(s));
                }

                ArrayList<Double>  sumSScoreSet = new ArrayList<>();
                sumSScore = sumSScore.replace("[", "");
                sumSScore = sumSScore.replace("]", "");
                for (String s :  sumSScore.split(",\\s") ){
                    sumSScoreSet.add(Double.valueOf(s));
                }

                numUIDSMappedQPhos.put(Long.valueOf(nbhdBaseDBID), numUIDsMappedSet);
                numEntitiesMappedQPhos.put(Long.valueOf(nbhdBaseDBID), nbhdMeasuredSet);
                avgSScoreQPhos.put(Long.valueOf(nbhdBaseDBID), avgSScoreSet);
                sumSScoreQPhos.put(Long.valueOf(nbhdBaseDBID), sumSScoreSet);
            }
        }


        ////////////////////////////////////////////// Generatign stats for this Mapped data //////////////////////////////////////////////
        try (Transaction tx = graphDb.beginTx()) {

            // get all labels and find one that matches the experiment
            Boolean experimentLabel = true;
            String supportScoreString = "";
            HashSet<String> experiments = new HashSet<>();
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            for (String property: allPropertyKeys){
                if(property.replace("SUPPORT_SCORE_", "").equalsIgnoreCase(experiment)){ // if experiment name is in properties
                    experimentLabel = false; // dont throw error
                    supportScoreString = property;
                }else if (property.contains("SUPPORT_SCORE_")){ // gather experiment names in db
                    String support_score_ = property.replace("SUPPORT_SCORE_", "");
                    experiments.add(support_score_);
                }
            }
            // throw exception if experiment name given is not in the database
            if(experimentLabel){
                throw new InputException(experiment + " is not currently in this database: " + databaseDir +
                        "\nExperiments in this database are: " + experiments);
            }

            // if we get this far the experiment is in the db
            String scoredByString = "SCORED_BY_" + experiment;
            String mappedString = "MAPPED_" + experiment;


            FileWriter fstream0 = new FileWriter(outputFile + "/NeighbourhoodAnalysis_" + experiment + ".tsv");
            BufferedWriter out0 = new BufferedWriter(fstream0);
            out0.write("nbhdBaseDBID\tnbhdBaseUID\tDisplayName\tNumIntegrated\tnbhdSize\tentitiesMapped\tentitiesMappedFDR\tnumUIDs\tnumUIDsMapped\tnumUIDsFDR\tavgSScore\tavgSScoreFDR\tsumSScore\tsumSScoreFDR\n");
            out0.close();
            //Create File to write inputs/Outputs to
            FileWriter fstream = new FileWriter(outputFile + "/NeighbourhoodAnalysis_" + experiment + ".tsv", true);
            BufferedWriter out = new BufferedWriter(fstream);

            // generate neighbourhoods and count stuff
            Integer countEntSig = 0;
            Integer countUIDSig = 0;
            Integer countAvgSig = 0;
            Integer countSumSig = 0;
            Integer countNumNbhds = 0;
            ResourceIterator<Node> allPeNodes = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            for (ResourceIterator<Node> it = allPeNodes; it.hasNext(); ) {
                Node physicalEntity = it.next();
                // if pe is measured - gen nbhd
                ArrayList<Double> nbhdSScores = new ArrayList<>();
                Integer nbhdSizeInt = 0;
                Integer nbhdMeasured = 0;
                Integer numIntegrated = 0;
                HashSet<String> nbhdUids = new HashSet<>();
                HashSet<String> numUIDsMeasured = new HashSet<>();
                if(physicalEntity.hasLabel(Label.label("Protein")) | physicalEntity.hasLabel(Label.label("Complex"))) {
                    HashSet<Node> nbhd = neighbourTraversal(physicalEntity, depth, graphDb);
                    countNumNbhds ++;


                    // count the total neighbourhood size
                    // count number of measured things
                    for (Node node : nbhd) {
                        if (node.hasLabel(Label.label("Protein")) | node.hasLabel(Label.label("Complex"))) {
                            nbhdSizeInt++;
                            //count integrated
                            if (node.hasProperty(PropertyType.INTEGRATED.toString())) {
                                numIntegrated++;
                            }
                        }

                        // get UID's from proteins and complexes in the network & count the # measured
                        if(node.hasLabel(Label.label("Protein"))){
                            Iterable<Relationship> uidRels = node.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                            for(Relationship uidRel: uidRels){
                                nbhdUids.add(uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                                if(uidRel.getStartNode().hasProperty(mappedString)){
                                    numUIDsMeasured.add(uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                                }
                            }
                        }
                        if( node.hasLabel(Label.label("Complex"))){
                            String cplxUIDs = node.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                            String[] split = cplxUIDs.split(",\\s");
                            for (String uid: split) {
                                Pattern p = Pattern.compile(UID_PATTERN);
                                Matcher m = p.matcher(uid);
                                if (m.find()) {
                                    nbhdUids.add(m.group(0));
                                }
                            }
                            if(node.hasProperty(mappedString)){
                                String mappedUIDs = node.getProperty(mappedString).toString();
                                String[] split1 = mappedUIDs.split(",\\s");
                                for (String uid: split1){
                                    numUIDsMeasured.add(uid);
                                }
                            }
                        }

                        // count measured
                        if (node.hasProperty(supportScoreString)) {
                            if (Double.parseDouble(node.getProperty(supportScoreString).toString()) > 0) {
                                nbhdMeasured++;
                            }
                            nbhdSScores.add(Double.parseDouble(node.getProperty(supportScoreString).toString()));
                        }
                    }

                    // AVG & SUM & STDDEV
                    Double ssAvg = 0.0;
                    Double ssSum = 0.0;
                    for (Integer i = 0; i < nbhdSScores.size(); i++) {
                        if (nbhdSScores.get(i) >= 0) {
                            ssAvg += nbhdSScores.get(i);
                            ssSum += nbhdSScores.get(i);
                        }
                    }
                    if (!ssSum.equals(0.0)) {
                        ssAvg = ssAvg / nbhdSScores.size();
                    }



                    ////////////// gathering outputs
                    String uid = "";

                    if (physicalEntity.hasProperty(scoredByString)) {
                        uid = physicalEntity.getProperty(scoredByString).toString();
                    } else if (physicalEntity.hasProperty(PropertyType.UNIPROT_ID.toString())) {
                        uid = physicalEntity.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                    } else {
                        Iterable<Relationship> uidRels = physicalEntity.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                        if (getLength(uidRels) == 1) {
                            for (Relationship uidRel : uidRels) {
                                uid = uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                            }
                        } else {
                            for (Relationship uidRel : uidRels) {
                                uid = uid + uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString() + ", ";
                            }
                        }
                    }

                    ////////////////////////////////////////////// Comparing stats for this Mapped data //////////////////////////////////////////////

                    Integer numUIDsMeasuredInt = numUIDsMeasured.size();
                    System.out.println(physicalEntity.getId());
                    ArrayList<Integer> integers = numUIDSMappedQPhos.get(physicalEntity.getId());
                    Integer countNumSmaller = 0;
                    for (Integer i : integers) {
                        if(numUIDsMeasuredInt <= i){ // actually counting things that are larger
                            countNumSmaller++;
                        }
                    }
                    Double numUIDsMeasuredFDR =  (countNumSmaller*1.0)/(integers.size()*1.0);

                    ArrayList<Integer> integers1 = numEntitiesMappedQPhos.get(physicalEntity.getId());
                    System.out.println(integers1);
                    System.out.println(nbhdMeasured);
                    Integer countNumSmaller1 = 0;
                    for (Integer i : integers1) {
                        if(nbhdMeasured <= i){
                            countNumSmaller1++;
                        }
                    }
                    Double entitesMeasuredFDR =  (countNumSmaller1*1.0)/(integers1.size()*1.0);
                    System.out.println(entitesMeasuredFDR);

                    ArrayList<Double> integers2 = avgSScoreQPhos.get(physicalEntity.getId());
                    Integer countNumSmaller2 = 0;
                    for (Double i : integers2) {
                        if(ssAvg <= i){
                            countNumSmaller2++;
                        }
                    }
                    Double ssAvgFDR =  (countNumSmaller2*1.0)/(integers2.size()*1.0);

                    ArrayList<Double> integers3 = sumSScoreQPhos.get(physicalEntity.getId());
                    Integer countNumSmaller3 = 0;
                    for (Double i : integers3) {
                        if(ssSum <= i){
                            countNumSmaller3++;
                        }
                    }
                    Double ssSumFDR =  (countNumSmaller3*1.0)/(integers3.size()*1.0);

                    if(entitesMeasuredFDR < 0.05){
                        countEntSig ++;
                    }
                    if(numUIDsMeasuredFDR < 0.05){
                        countUIDSig ++;
                    }
                    if(ssAvgFDR < 0.05){
                        countAvgSig ++;
                    }
                    if(ssSumFDR < 0.05){
                        countSumSig ++;
                    }

                    out.write(physicalEntity.getId() + "\t"
                            + uid
                            + "\t" + physicalEntity.getProperty(PropertyType.DISPLAY_NAME.toString())
                            + "\t" + numIntegrated
                            + "\t" + nbhdSizeInt
                            + "\t" + nbhdMeasured
                            + "\t" +entitesMeasuredFDR
                            + "\t" + nbhdUids.size()
                            + "\t" + numUIDsMeasured.size()
                            + "\t" +numUIDsMeasuredFDR
                            + "\t " + ssAvg
                            + "\t" +ssAvgFDR
                            + "\t" + ssSum
                            + "\t" +ssSumFDR
                            + "\n");

                    // this is where you add nbhd stats to nodes
                    physicalEntity.setProperty(("NBHD_entity_measured_FDR_" + experiment), entitesMeasuredFDR);
                    physicalEntity.setProperty(("NBHD_uids_measured_FDR_" + experiment), numUIDsMeasuredFDR);
                    physicalEntity.setProperty(("NBHD_avg_FDR_" + experiment), ssAvgFDR);
                    physicalEntity.setProperty(("NBHD_sum_FDR_" + experiment), ssSumFDR);
                }
            }
            System.out.println("Num Sig nbhds w Entities: " + countEntSig);
            System.out.println("Num Sig nbhds w UIDs: " + countUIDSig);
            System.out.println("Num Sig nbhds w AVG: " + countAvgSig);
            System.out.println("Num Sig nbhds w SUM: " + countSumSig);
            System.out.println("Total number of nbhds: " + countNumNbhds);
            tx.success();
            out.close();
        }

    }

    void nbhdAnalysis(Integer depth, String experiment, File empiricalDistribution1, File empiricalDistribution2) throws IOException {

        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        ////////////////////////////////////////////// Reading in empirical dists //////////////////////////////////////////////
        BufferedReader BR1 = new BufferedReader(new FileReader(empiricalDistribution1));
        BufferedReader BR2 = new BufferedReader(new FileReader(empiricalDistribution1));
        String line = "";

        Integer nbhdBaseDBIDCol =0;
        Integer numEntitiesCol = 0;
        Integer numIntegratedCol = 0;
        Integer numUIDsCol = 0;
        Integer numUIDsMappedCol = 0;
        Integer nbhdMesuredCol = 0;
        Integer avgSScoreCol = 0;
        Integer sumSScoreCol = 0;

        System.out.println("Reading in all distributions per Neighbourhood");

        Integer countID = 0;

        HashMap<Long, ArrayList<Integer>> numUIDSMappedQPhos = new HashMap<>();
        HashMap<Long, ArrayList<Integer>> numEntitiesMappedQPhos = new HashMap<>();
        HashMap<Long, ArrayList<Double>> avgSScoreQPhos = new HashMap<>();
        HashMap<Long, ArrayList<Double>> sumSScoreQPhos = new HashMap<>();

        while ((line = BR1.readLine()) != null) {
            String[] columns = line.split("\t");
            Integer colLen = columns.length;
            if (countID == 0) { // if first line in file
                for (int i = 0; i < colLen; i++) {
                    if (columns[i].equalsIgnoreCase("nbhdBaseDBID")) {
                        nbhdBaseDBIDCol = i;
                    }
                    if (columns[i].equals("numEntities")) {
                        numEntitiesCol = i;
                    }
                    if (columns[i].equals("NumIntegrated")) {
                        numIntegratedCol = i;
                    }
                    if (columns[i].equals("numUIDs")) {
                        numUIDsCol = i;
                    }
                    if (columns[i].equals("numUIDsMapped")) {
                        numUIDsMappedCol = i;
                    }
                    if (columns[i].equals("nbhdMeasured")) {
                        nbhdMesuredCol = i;
                    }
                    if (columns[i].equals("avgSScore")) {
                        avgSScoreCol = i;
                    }
                    if (columns[i].equals("sumSScore")) {
                        sumSScoreCol = i;
                    }

                }
                countID++;
            } else { //for all other lines

                String nbhdBaseDBID = columns[nbhdBaseDBIDCol];
                String numEntities = columns[numEntitiesCol];
                String numIntegrated = columns[numIntegratedCol];
                String numUIDs = columns[numUIDsCol];
                String numUIDsMapped = columns[numUIDsMappedCol];
                String nbhdMesured = columns[nbhdMesuredCol];
                String avgSScore = columns[avgSScoreCol];
                String sumSScore = columns[sumSScoreCol];

                numUIDsMapped = numUIDsMapped.replace("[", "");
                numUIDsMapped = numUIDsMapped.replace("]", "");
                ArrayList<Integer>  numUIDsMappedSet = new ArrayList<>();
                for (String s :  numUIDsMapped.split(",\\s") ){
                    numUIDsMappedSet.add(Integer.valueOf(s));
                }

                ArrayList<Integer>  nbhdMeasuredSet = new ArrayList<>();
                nbhdMesured = nbhdMesured.replace("[", "");
                nbhdMesured = nbhdMesured.replace("]", "");
                for (String s :  nbhdMesured.split(",\\s") ){
                    nbhdMeasuredSet.add(Integer.valueOf(s));
                }

                ArrayList<Double>  avgSScoreSet = new ArrayList<>();
                avgSScore = avgSScore.replace("[", "");
                avgSScore = avgSScore.replace("]", "");
                for (String s :  avgSScore.split(",\\s") ){
                    avgSScoreSet.add(Double.valueOf(s));
                }

                ArrayList<Double>  sumSScoreSet = new ArrayList<>();
                sumSScore = sumSScore.replace("[", "");
                sumSScore = sumSScore.replace("]", "");
                for (String s :  sumSScore.split(",\\s") ){
                    sumSScoreSet.add(Double.valueOf(s));
                }

                numUIDSMappedQPhos.put(Long.valueOf(nbhdBaseDBID), numUIDsMappedSet);
                numEntitiesMappedQPhos.put(Long.valueOf(nbhdBaseDBID), nbhdMeasuredSet);
                avgSScoreQPhos.put(Long.valueOf(nbhdBaseDBID), avgSScoreSet);
                sumSScoreQPhos.put(Long.valueOf(nbhdBaseDBID), sumSScoreSet);
            }
        }

        line = "";

        nbhdBaseDBIDCol =0;
        numEntitiesCol = 0;
        numIntegratedCol = 0;
        numUIDsCol = 0;
        numUIDsMappedCol = 0;
        nbhdMesuredCol = 0;
        avgSScoreCol = 0;
        sumSScoreCol = 0;

        countID = 0;
        while ((line = BR2.readLine()) != null) {
            String[] columns = line.split("\t");
            Integer colLen = columns.length;
            if (countID == 0) { // if first line in file
                for (int i = 0; i < colLen; i++) {
                    if (columns[i].equalsIgnoreCase("nbhdBaseDBID")) {
                        nbhdBaseDBIDCol = i;
                    }
                    if (columns[i].equals("numEntities")) {
                        numEntitiesCol = i;
                    }
                    if (columns[i].equals("NumIntegrated")) {
                        numIntegratedCol = i;
                    }
                    if (columns[i].equals("numUIDs")) {
                        numUIDsCol = i;
                    }
                    if (columns[i].equals("numUIDsMapped")) {
                        numUIDsMappedCol = i;
                    }
                    if (columns[i].equals("nbhdMeasured")) {
                        nbhdMesuredCol = i;
                    }
                    if (columns[i].equals("avgSScore")) {
                        avgSScoreCol = i;
                    }
                    if (columns[i].equals("sumSScore")) {
                        sumSScoreCol = i;
                    }

                }
                countID++;
            } else { //for all other lines

                String nbhdBaseDBID = columns[nbhdBaseDBIDCol];
                String numEntities = columns[numEntitiesCol];
                String numIntegrated = columns[numIntegratedCol];
                String numUIDs = columns[numUIDsCol];
                String numUIDsMapped = columns[numUIDsMappedCol];
                String nbhdMesured = columns[nbhdMesuredCol];
                String avgSScore = columns[avgSScoreCol];
                String sumSScore = columns[sumSScoreCol];

                numUIDsMapped = numUIDsMapped.replace("[", "");
                numUIDsMapped = numUIDsMapped.replace("]", "");
                ArrayList<Integer>  numUIDsMappedSet = new ArrayList<>();
                for (String s :  numUIDsMapped.split(",\\s") ){
                    numUIDsMappedSet.add(Integer.valueOf(s));
                }

                ArrayList<Integer>  nbhdMeasuredSet = new ArrayList<>();
                nbhdMesured = nbhdMesured.replace("[", "");
                nbhdMesured = nbhdMesured.replace("]", "");
                for (String s :  nbhdMesured.split(",\\s") ){
                    nbhdMeasuredSet.add(Integer.valueOf(s));
                }

                ArrayList<Double>  avgSScoreSet = new ArrayList<>();
                avgSScore = avgSScore.replace("[", "");
                avgSScore = avgSScore.replace("]", "");
                for (String s :  avgSScore.split(",\\s") ){
                    avgSScoreSet.add(Double.valueOf(s));
                }

                ArrayList<Double>  sumSScoreSet = new ArrayList<>();
                sumSScore = sumSScore.replace("[", "");
                sumSScore = sumSScore.replace("]", "");
                for (String s :  sumSScore.split(",\\s") ){
                    sumSScoreSet.add(Double.valueOf(s));
                }

                // this time add to existing maps
                if(numUIDSMappedQPhos.containsKey(Long.valueOf(nbhdBaseDBID))){
                    ArrayList<Integer> integers = numUIDSMappedQPhos.get(Long.valueOf(nbhdBaseDBID));
                    integers.addAll(numUIDsMappedSet);
                    numUIDSMappedQPhos.put(Long.valueOf(nbhdBaseDBID), integers);
                }

                if(numEntitiesMappedQPhos.containsKey(Long.valueOf(nbhdBaseDBID))){
                    ArrayList<Integer> integers = numEntitiesMappedQPhos.get(Long.valueOf(nbhdBaseDBID));
                    integers.addAll(nbhdMeasuredSet);
                    numEntitiesMappedQPhos.put(Long.valueOf(nbhdBaseDBID), integers);
                }

                if(avgSScoreQPhos.containsKey(Long.valueOf(nbhdBaseDBID))){
                    ArrayList<Double> doubles = avgSScoreQPhos.get(Long.valueOf(nbhdBaseDBID));
                    doubles.addAll(avgSScoreSet);
                    avgSScoreQPhos.put(Long.valueOf(nbhdBaseDBID), doubles);
                }

                if(sumSScoreQPhos.containsKey(Long.valueOf(nbhdBaseDBID))){
                    ArrayList<Double> doubles = sumSScoreQPhos.get(Long.valueOf(nbhdBaseDBID));
                    doubles.addAll(sumSScoreSet);
                    sumSScoreQPhos.put(Long.valueOf(nbhdBaseDBID), doubles);
                }

            }
        }

        System.out.println(numUIDSMappedQPhos);
        System.out.println(numEntitiesMappedQPhos);
        System.out.println(avgSScoreQPhos);
        System.out.println(sumSScoreQPhos);


        ////////////////////////////////////////////// Generatign stats for this Mapped data //////////////////////////////////////////////
        try (Transaction tx = graphDb.beginTx()) {

            // get all labels and find one that matches the experiment
            Boolean experimentLabel = true;
            String supportScoreString = "";
            HashSet<String> experiments = new HashSet<>();
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            for (String property: allPropertyKeys){
                if(property.replace("SUPPORT_SCORE_", "").equalsIgnoreCase(experiment)){ // if experiment name is in properties
                    experimentLabel = false; // dont throw error
                    supportScoreString = property;
                }else if (property.contains("SUPPORT_SCORE_")){ // gather experiment names in db
                    String support_score_ = property.replace("SUPPORT_SCORE_", "");
                    experiments.add(support_score_);
                }
            }
            // throw exception if experiment name given is not in the database
            if(experimentLabel){
                throw new InputException(experiment + " is not currently in this database: " + databaseDir +
                        "\nExperiments in this database are: " + experiments);
            }

            // if we get this far the experiment is in the db
            String scoredByString = "SCORED_BY_" + experiment;
            String mappedString = "MAPPED_" + experiment;


            FileWriter fstream0 = new FileWriter(outputFile + "/NeighbourhoodAnalysis_" + experiment + ".tsv");
            BufferedWriter out0 = new BufferedWriter(fstream0);
            out0.write("nbhdBaseDBID\tnbhdBaseUID\tDisplayName\tNumIntegrated\tnbhdSize\tentitiesMapped\tentitiesMappedFDR\tnumUIDs\tnumUIDsMapped\tnumUIDsFDR\tavgSScore\tavgSScoreFDR\tsumSScore\tsumSScoreFDR\n");
            out0.close();
            //Create File to write inputs/Outputs to
            FileWriter fstream = new FileWriter(outputFile + "/NeighbourhoodAnalysis_" + experiment + ".tsv", true);
            BufferedWriter out = new BufferedWriter(fstream);

            // generate neighbourhoods and count stuffJ.
            Integer countEntSig = 0;
            Integer countUIDSig = 0;
            Integer countAvgSig = 0;
            Integer countSumSig = 0;
            Integer countNumNbhds = 0;
            ResourceIterator<Node> allPeNodes = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            for (ResourceIterator<Node> it = allPeNodes; it.hasNext(); ) {
                Node physicalEntity = it.next();
                // if pe is measured - gen nbhd
                ArrayList<Double> nbhdSScores = new ArrayList<>();
                Integer nbhdSizeInt = 0;
                Integer nbhdMeasured = 0;
                Integer numIntegrated = 0;
                HashSet<String> nbhdUids = new HashSet<>();
                HashSet<String> numUIDsMeasured = new HashSet<>();
                if(physicalEntity.hasLabel(Label.label("Protein")) | physicalEntity.hasLabel(Label.label("Complex"))) {
                    HashSet<Node> nbhd = neighbourTraversal(physicalEntity, depth, graphDb);
                    countNumNbhds ++;


                    // count the total neighbourhood size
                    // count number of measured things
                    for (Node node : nbhd) {
                        if (node.hasLabel(Label.label("Protein")) | node.hasLabel(Label.label("Complex"))) {
                            nbhdSizeInt++;
                            //count integrated
                            if (node.hasProperty(PropertyType.INTEGRATED.toString())) {
                                numIntegrated++;
                            }
                        }

                        // get UID's from proteins and complexes in the network & count the # measured
                        if(node.hasLabel(Label.label("Protein"))){
                            Iterable<Relationship> uidRels = node.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                            for(Relationship uidRel: uidRels){
                                nbhdUids.add(uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                                if(uidRel.getStartNode().hasProperty(mappedString)){
                                    numUIDsMeasured.add(uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString());
                                }
                            }
                        }
                        if( node.hasLabel(Label.label("Complex"))){
                            String cplxUIDs = node.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                            String[] split = cplxUIDs.split(",\\s");
                            for (String uid: split) {
                                Pattern p = Pattern.compile(UID_PATTERN);
                                Matcher m = p.matcher(uid);
                                if (m.find()) {
                                    nbhdUids.add(m.group(0));
                                }
                            }
                            if(node.hasProperty(mappedString)){
                                String mappedUIDs = node.getProperty(mappedString).toString();
                                String[] split1 = mappedUIDs.split(",\\s");
                                for (String uid: split1){
                                    numUIDsMeasured.add(uid);
                                }
                            }
                        }

                        // count measured
                        if (node.hasProperty(supportScoreString)) {
                            if (Double.parseDouble(node.getProperty(supportScoreString).toString()) > 0) {
                                nbhdMeasured++;
                            }
                            nbhdSScores.add(Double.parseDouble(node.getProperty(supportScoreString).toString()));
                        }
                    }

                    // AVG & SUM & STDDEV
                    Double ssAvg = 0.0;
                    Double ssSum = 0.0;
                    for (Integer i = 0; i < nbhdSScores.size(); i++) {
                        if (nbhdSScores.get(i) >= 0) {
                            ssAvg += nbhdSScores.get(i);
                            ssSum += nbhdSScores.get(i);
                        }
                    }
                    if (!ssSum.equals(0.0)) {
                        ssAvg = ssAvg / nbhdSScores.size();
                    }



                    ////////////// gathering outputs
                    String uid = "";

                    if (physicalEntity.hasProperty(scoredByString)) {
                        uid = physicalEntity.getProperty(scoredByString).toString();
                    } else if (physicalEntity.hasProperty(PropertyType.UNIPROT_ID.toString())) {
                        uid = physicalEntity.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                    } else {
                        Iterable<Relationship> uidRels = physicalEntity.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                        if (getLength(uidRels) == 1) {
                            for (Relationship uidRel : uidRels) {
                                uid = uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                            }
                        } else {
                            for (Relationship uidRel : uidRels) {
                                uid = uid + uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString() + ", ";
                            }
                        }
                    }

                    ////////////////////////////////////////////// Comparing stats for this Mapped data //////////////////////////////////////////////

                    Integer numUIDsMeasuredInt = numUIDsMeasured.size();
                    ArrayList<Integer> integers = numUIDSMappedQPhos.get(physicalEntity.getId());
                    Integer countNumSmaller = 0;
                    for (Integer i : integers) {
                        if(numUIDsMeasuredInt <= i){
                            countNumSmaller++;
                        }
                    }
                    Double numUIDsMeasuredFDR =  (countNumSmaller*1.0)/(integers.size()*1.0);

                    ArrayList<Integer> integers1 = numEntitiesMappedQPhos.get(physicalEntity.getId());
                    Integer countNumSmaller1 = 0;
                    for (Integer i : integers1) {
                        if(nbhdMeasured <= i){
                            countNumSmaller1++;
                        }
                    }
                    Double entitesMeasuredFDR =  (countNumSmaller1*1.0)/(integers1.size()*1.0);

                    ArrayList<Double> integers2 = avgSScoreQPhos.get(physicalEntity.getId());
                    Integer countNumSmaller2 = 0;
                    for (Double i : integers2) {
                        if(ssAvg <= i){
                            countNumSmaller2++;
                        }
                    }
                    Double ssAvgFDR =  (countNumSmaller2*1.0)/(integers2.size()*1.0);

                    ArrayList<Double> integers3 = sumSScoreQPhos.get(physicalEntity.getId());
                    Integer countNumSmaller3 = 0;
                    for (Double i : integers3) {
                        if(ssSum <= i){
                            countNumSmaller3++;
                        }
                    }
                    Double ssSumFDR =  (countNumSmaller3*1.0)/(integers3.size()*1.0);

                    if(entitesMeasuredFDR < 0.05){
                        countEntSig ++;
                    }
                    if(numUIDsMeasuredFDR < 0.05){
                        countUIDSig ++;
                    }
                    if(ssAvgFDR < 0.05){
                        countAvgSig ++;
                    }
                    if(ssSumFDR < 0.05){
                        countSumSig ++;
                    }

                    out.write(physicalEntity.getId() + "\t"
                            + uid
                            + "\t" + physicalEntity.getProperty(PropertyType.DISPLAY_NAME.toString())
                            + "\t" + numIntegrated
                            + "\t" + nbhdSizeInt
                            + "\t" + nbhdMeasured
                            + "\t" +entitesMeasuredFDR
                            + "\t" + nbhdUids.size()
                            + "\t" + numUIDsMeasured.size()
                            + "\t" +numUIDsMeasuredFDR
                            + "\t " + ssAvg
                            + "\t" +ssAvgFDR
                            + "\t" + ssSum
                            + "\t" +ssSumFDR
                            + "\n");

                    // this is where you add nbhd stats to nodes
                    physicalEntity.setProperty(("NBHD_entity_measured_FDR_" + experiment), entitesMeasuredFDR);
                    physicalEntity.setProperty(("NBHD_uids_measured_FDR_" + experiment), numUIDsMeasuredFDR);
                    physicalEntity.setProperty(("NBHD_avg_FDR_" + experiment), ssAvgFDR);
                    physicalEntity.setProperty(("NBHD_sum_FDR_" + experiment), ssSumFDR);
                }

            }
            System.out.println("Num Sig nbhds w Entities: " + countEntSig);
            System.out.println("Num Sig nbhds w UIDs: " + countUIDSig);
            System.out.println("Num Sig nbhds w AVG: " + countAvgSig);
            System.out.println("Num Sig nbhds w SUM: " + countSumSig);
            System.out.println("Total number of nbhds: " + countNumNbhds);
            tx.success();
            out.close();
        }

    }

    /**
     * A helper function that takes a peptide, it's id, and the ID2Start hashmap and
     * returns a set of phosphorylations on that peptide
     * @param ID
     * @param pep
     * @param ID2Start
     * @return
     */
    private HashSet<String> generateModList(int ID, String pep, HashMap<String, String> ID2Start){

        //if in this format -rm excess _(ac)AAAITDM(ox)ADLEELSRLS(ph)PLPPGS(ph)PGSAAR_
        Matcher m = Pattern.compile("\\(..\\)").matcher(pep);
        if (m.find()){
            if (m.group(0).equals("(ac)") | m.group(0).equals("(ox)") | m.group(0).equals("(ph)") ){
                pep = pep.replaceAll("\\(ac\\)|\\(ox\\)|\\_", "");
            }
        }

        HashSet<String> mods = new HashSet<>();
        int pcount = 0; //# phopshorylations in this pep

        //if in this format _AAAITDMADLEELSRLS(ph)PLPPGS(ph)PGSAAR_
        if(pep.contains("(ph)")){
            for(int i = 0; i <pep.length();i++ ){
                if(pep.charAt(i) == '('){
                    pcount ++;
                    char modType = pep.charAt((i-1)); // bc the char before is the phosphorylated one
                    String x = ID2Start.get(String.valueOf(ID));
                    int modLoc = 0;
                    if(pcount == 1){
                        modLoc = (i-1) + Integer.valueOf(x); // if no '(ph)' before just its the last aa + start of pep
                    }else{
                        modLoc = ((i-1) + Integer.valueOf(x)) - ((pcount -1)*4); // if multi '(ph)'s subtract from pep length
                    }
                    String modTypeStr = String.valueOf(modType);
                    String modLocStr = String.valueOf(modLoc);
                    String mod = modTypeStr +"_"+ modLocStr;
                    mods.add(mod);
                }
            }
        }else{ //if in this format AAAITDMADLEELSRLpSPLPPGpSPGSAAR
            for(int i = 0; i <pep.length();i++ ){
                if (pep.charAt(i) == 'p'){
                    pcount ++;
                    char modType = pep.charAt(i+1); // get character after
                    String x = ID2Start.get(String.valueOf(ID)); // get start of pep
                    int modLoc = (i-pcount) + Integer.valueOf(x) + 1; // subtract the number of pho's already passed from current count and add the start of the pep + 1 to account for counting from 0
                    String modTypeStr = String.valueOf(modType);
                    String modLocStr = String.valueOf(modLoc);
                    String mod = modTypeStr +"_"+ modLocStr;
                    mods.add(mod);
                }
            }
        }

        return mods;
    }

    /**
     * a helper function that takes in a fasta file, a UID and the modfied peptide
     * returns the amino acid the peptide starts on
     * @param fasta
     * @param LRP
     * @param modPep
     * @return
     */
    private String generateStartOfPep(HashMap<String, String> fasta, String LRP, String modPep){
        String start = "";

        // get fasta seq
        String seq = "";
        if(fasta.containsKey(LRP)){
            seq = fasta.get(LRP);
        }

        // get mod pep
        modPep = modPep.replaceAll("\\(..\\)|[p]|\\_", "");

        int i= 0;
        if(seq.contains(modPep)){
            // find in seq
            i = seq.indexOf(modPep);
            i = i+1;
            start = String.valueOf(i);
        }else{
            //System.out.println(LRP +"'s fasta doesnt contain "+ modPep);
        }

        return start;

    }

    /**
     * a helper function that takes in the starting point of the peptide, and the modified peptide
     * returns the amino acid the peptide ends on
     * @param start
     * @param modPep
     * @return
     */
    private String generateEndOfPep(String start, String modPep){
        String end = "";

        if(!start.isEmpty()){
            int strt = Integer.valueOf(start);
            modPep = modPep.replaceAll("\\(..\\)|[p]|\\_", "");
            int length = modPep.length();
            int nd = strt + length-1;
            end = String.valueOf(nd);
        }

        return end;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // getNeighbours gets neighbourhoods for specific porteins or files of proteins

    /**
     * overloaded - Takes in a file of uids or db ids and passes them to getNeighbours
     * @param pathToProt
     * @param depth
     * @throws IOException
     */
    public void getNeighbours(File pathToProt, Integer depth) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        File OUTPUT_PATH = outputFile;
        File NEIGHBOURHOOD_PATH = new File(OUTPUT_PATH + "/neighbourhoods");

        NEIGHBOURHOOD_PATH.mkdirs();

        // get Proteins from File
        BufferedReader BR = null;
        String line = "";
        String csvSplitBy = ",";
        BR = new BufferedReader(new FileReader(pathToProt));
        Set<String> UIDList = new HashSet<>();
        while ((line = BR.readLine()) != null) {
            UIDList.add(line);
        }
        for (String uid: UIDList){

            getNeighbours(uid,depth);
        }

    }
    /**
     * Takes in a sting and determines if it's a uid or db_id and then get it's neighbourhood
     * @param uid
     * @param depth
     */
    public void getNeighbours(String uid, Integer depth) {
        System.out.println("Depth "+ depth);
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();

        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        File OUTPUT_PATH = outputFile;

        try (Transaction tx = graphDb.beginTx()) {


            Pattern p = Pattern.compile(UID_PATTERN);
            Matcher m = p.matcher(uid);
            if (m.find()) {     //////////////////////////// if its a uniprotid
                String theGroup = m.group(0);
                Node node = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()),
                        PropertyType.UNIPROT_ID.toString(), theGroup);
                Iterable<Relationship> relationships = node.getRelationships(RelTypes.ID_BELONGS_TO);
                for (Relationship relationship : relationships) {
                    Node prot = relationship.getEndNode();
                    neighbourTraversal(node,depth,graphDb);

                }
            }else{              //////////////////////////// its a node id
                Node node = graphDb.getNodeById(Long.valueOf(uid));
                if (node == null){
                    throw new Exception("Node "+ uid+" not found");
                }else {
                    neighbourTraversal(node, depth, graphDb);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        graphDb.shutdown();
    }

    ////////////////////////////

    // gets inputs and outputs of network
    private void getIOs() throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        File OUTPUT_PATH = outputFile;

        try (Transaction tx = graphDb.beginTx()) {

            //Create File to write inputs/Outputs to
            FileWriter fstream;
            BufferedWriter out;
            // create your filewriter and bufferedreader
            fstream = new FileWriter(OUTPUT_PATH + "/IOs.tsv");
            out = new BufferedWriter(fstream);

            out.write("ID\tIO\n");

            // get all nodes tagged physical entity
            ResourceIterator<Node> physicalEntity = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            while (physicalEntity.hasNext()) {

                Node pe = physicalEntity.next();
                Iterable<Relationship> OUTrelationships = pe.getRelationships(RelTypes.OUTPUT);
                Integer outLen = getLength(OUTrelationships);

                Iterable<Relationship> INrelationships = pe.getRelationships(RelTypes.INPUT);
                Iterable<Relationship> CTLrelationships = pe.getRelationships(RelTypes.CONTROLS);
                Integer inLen = getLength(INrelationships) + getLength(CTLrelationships);

                if (inLen == 0) {
                    pe.addLabel(Label.label("OUTPUT"));
                    out.write(pe.getId() + "\t" + "OUTPUT" + "\n");
                }
                if (outLen == 0) {
                    pe.addLabel(Label.label("INPUT"));
                    out.write(pe.getId() + "\t" + "INPUT" + "\n");
                }

            }

            tx.success();

            out.close();
        }
        graphDb.shutdown();

    }

    ////////////////////////////

    // gets inputs and outputs of a protein of interest

    void getXPuts(File pathToUIDs, String io_decision_string) throws IOException {
        // get Proteins from File
        BufferedReader BR = null;
        String line = "";
        String csvSplitBy = ",";
        BR = new BufferedReader(new FileReader(pathToUIDs));
        Set<String> UIDList = new HashSet<>();
        while ((line = BR.readLine()) != null) {
            UIDList.add(line);
        }

        for(String uid :UIDList){
            getXPuts(uid ,io_decision_string );
        }

    }
    void getXPuts(String uid, String io_decision_string) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        File OUTPUT_PATH = outputFile;



        //Create File to write inputs to
        FileWriter fstream;
        BufferedWriter out;
        // create your filewriter and bufferedreader
        fstream = new FileWriter(OUTPUT_PATH + "/IOtoNodes.tsv");
        out = new BufferedWriter(fstream);


        io_decision_string = io_decision_string.toLowerCase();

        if (io_decision_string.equalsIgnoreCase("inputs") | io_decision_string.equalsIgnoreCase("i")) {
            out.write("ID\tUID\tRXMID\tInputID\tInputRXMID\n");
            try (Transaction tx = graphDb.beginTx()) {
                // get all of its uid nodes
                Node node = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), uid);
                // if it doesnt exist continue
                if (node == null) {
                    node = graphDb.getNodeById(Long.valueOf(uid));
                    if(node == null){
                        throw new Exception("Node "+ uid+" not found");
                    }
                }
                // node not found
                // get it's proteins
                Iterable<Relationship> relationships = node.getRelationships(RelTypes.ID_BELONGS_TO);
                // for each protein it maps to
                for (Relationship relationship : relationships) {
                    Long end = relationship.getEndNodeId();
                    Node prot = graphDb.getNodeById(end);
                    Iterable<Relationship> PROTrelationships = prot.getRelationships(RelTypes.OUTPUT);
                    for (Relationship PROTrelationship : PROTrelationships) {
                        out.write(prot.getId() + "\t" +
                                uid + "\t" +
                                prot.getProperty(PropertyType.DISPLAY_NAME.toString()) + "\t" +
                                PROTrelationship.getStartNodeId() + "\t" +
                                PROTrelationship.getStartNode().getProperty(PropertyType.DISPLAY_NAME.toString()) + "\n");
                    }
                }


                tx.success();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        } else if (io_decision_string.equalsIgnoreCase("outputs") || io_decision_string.equalsIgnoreCase("o")) {
            out.write("ID\tUID\tRXMID\tOutputID\tOutputRXMID\n");
            try (Transaction tx = graphDb.beginTx()) {

                // get all of its uid nodes
                Node node = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), uid);
                // if it doesnt exist continue
                if (node == null) {
                    node = graphDb.getNodeById(Long.valueOf(uid));
                }else{
                    // node not found
                    throw new Exception("Node "+ uid+" not found");
                }
                // if it exists get all if its inputs
                // get it's proteins
                Iterable<Relationship> relationships = node.getRelationships(RelTypes.ID_BELONGS_TO);
                // for each protein it maps to
                for (Relationship relationship : relationships) {
                    Long end = relationship.getEndNodeId();
                    Node prot = graphDb.getNodeById(end);
                    Iterable<Relationship> PROTrelationships = prot.getRelationships(RelTypes.INPUT);
                    for (Relationship PROTrelationship : PROTrelationships) {
                        out.write(prot.getId() + "\t" +
                                uid + "\t" +
                                prot.getId() + "\t" +
                                PROTrelationship.getEndNodeId() + "\t" +
                                PROTrelationship.getEndNode().getId() + "\n");
                    }
                    Iterable<Relationship> PROTrelationships2 = prot.getRelationships(RelTypes.CONTROLS);
                    for (Relationship PROTrelationship : PROTrelationships2) {
                        out.write(prot.getId() + "\t" +
                                uid + "\t" +
                                prot.getId() + "\t" +
                                PROTrelationship.getEndNodeId() + "\t" +
                                PROTrelationship.getEndNode().getId()+ "\n");
                    }
                }

                tx.success();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        } else if (io_decision_string.equalsIgnoreCase("io") ) {
            out.write("ID\tUID\tRXMID\tType\tID\tRXMID\n");
            try (Transaction tx = graphDb.beginTx()) {

                // get all of its uid nodes
                Node node = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), uid);
                // if it doesnt exist continue
                if (node == null) {
                    node = graphDb.getNodeById(Long.valueOf(uid));
                }else{
                    // node not found
                    throw new Exception("Node "+ uid+" not found");
                }
                Iterable<Relationship> relationships = node.getRelationships(RelTypes.ID_BELONGS_TO);
                // for each protein it maps to
                for (Relationship relationship : relationships) {
                    Long end = relationship.getEndNodeId();
                    Node prot = graphDb.getNodeById(end);
                    Iterable<Relationship> PROTrelationships = prot.getRelationships();
                    for (Relationship PROTrelationship : PROTrelationships) {
                        if (PROTrelationship.isType(RelTypes.INPUT) || PROTrelationship.isType(RelTypes.CONTROLS)) {
                            out.write(prot.getId() + "\t" +
                                    uid + "\t" +
                                    prot.getId() +
                                    "\tOUTPUTS:\t" +
                                    PROTrelationship.getEndNode() + "\t" +
                                    PROTrelationship.getEndNode().getId() + "\n");

                        } else if (PROTrelationship.isType(RelTypes.OUTPUT)) {
                            out.write(prot.getId() + "\t" +
                                    uid + "\t" + prot.getId() +
                                    "\tINPUTS:\t" +
                                    PROTrelationship.getStartNode() + "\t" +
                                    PROTrelationship.getStartNode().getId() + "\n");
                        } else {
                            continue;
                        }
                    }
                }


                tx.success();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        } else {
            throw new IllegalArgumentException("io_decision_string must equal Input, Output or IO");
        }

        out.close();
    } // TODO FIX THIS ! ALL WRONG !
    ////////////////////////////

    // null dist for tool

    private static void registerShutdownHook( final GraphDatabaseService graphDB ) {
        // Registers a shutdown hook for the Neo4j instance so that it
        // shuts down nicely when the VM exits (even if you "Ctrl-C" the
        // running application).
        Runtime.getRuntime().addShutdownHook( new Thread()
        {
            @Override
            public void run()
            {
                graphDB.shutdown();
            }
        } );
    }


}
