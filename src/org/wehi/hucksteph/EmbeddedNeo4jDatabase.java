package org.wehi.hucksteph;

import org.biopax.paxtools.model.level2.complex;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.graphdb.traversal.*;
import org.neo4j.unsafe.impl.batchimport.input.InputException;
import scala.Int;


import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class EmbeddedNeo4jDatabase {

    private File databaseDir;
    private File outputFile;


    private final String UID_PATTERN = "[OPQ][0-9][A-Z0-9]{3}[0-9](\\-[0-9*]{1,2})?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(\\-[0-9*]{1,2})?";
    private final String humanUniProt = "https://www.uniprot.org/uniprot/?query=organism:9606&format=fasta&include=yes";
    private final String mouseUniProt = "https://www.uniprot.org/uniprot/?query=organism:10090&format=fasta&include=yes";

    public File getDatabaseDir() {
        return databaseDir;
    }

    public File getOutputFile() {
        return outputFile;
    }

    public EmbeddedNeo4jDatabase(File databaseDir, File outputFile) {
        this.databaseDir = databaseDir;
        this.outputFile = outputFile;

    }

    public EmbeddedNeo4jDatabase(File databaseDir) {
        this.databaseDir = databaseDir;
    }

    /**
     * Prints all nodes and all relationships in a database
     */
    public void printDatabase() {
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterable<Node> nodeIterable = graphDb.getAllNodes();
            System.out.println("\nNodes:\n");
            for (Node node : nodeIterable) {
                System.out.println(node.getAllProperties());
            }

            ResourceIterable<Relationship> allRelationships = graphDb.getAllRelationships();
            System.out.println("\nRelationships:\n");
            for(Relationship relationship: allRelationships){
                String startNode = "";
                String endNode = "";

                if(relationship.getStartNode().hasProperty(PropertyType.DISPLAY_NAME.toString())){
                    startNode = relationship.getStartNode().getProperty(PropertyType.DISPLAY_NAME.toString()).toString();
                }else{
                    startNode = relationship.getStartNode().getId() + "";
                }
                if(relationship.getStartNode().hasProperty(PropertyType.DISPLAY_NAME.toString())){
                    endNode = relationship.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).toString();
                }else{
                    endNode = relationship.getEndNode().getId() + "";
                }

                System.out.println(startNode +
                        " -> "+
                        relationship.getType().name()+
                        relationship.getAllProperties().toString()+
                        " -> " +
                        endNode);
            }
            tx.success();
        }
        graphDb.shutdown();
    }

    /**
     * Gets the species node and assigns the species to the global private string
     * @param graphDb
     */
    public String getSpecies(GraphDatabaseService graphDb){

        String species = "";
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> speciesNodes = graphDb.findNodes(Label.label("SPECIES"));

            while(speciesNodes.hasNext()){
                Node speciesNode = speciesNodes.next();
                String speciesStr = speciesNode.getProperty("Species").toString();
                if(speciesStr.equalsIgnoreCase("human")){
                    species = "Human";
                }else{
                    species = "Mouse";
                }
            }
            tx.success();
        }
        return species;
    }

    /**
     * prints all UniProt IDs in database including deleted IDs
     * @throws IOException
     */
    public void printAllUniProtIDs() throws IOException {
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        File uniprotFile = new File(outputFile + "/UniProtIDs.tsv");
        uniprotFile.createNewFile();
        if(!uniprotFile.createNewFile()){
            //TODO if file exists ask if you can write over
            System.out.println("UniProtIDs.tsv File already exists\nWriting over ");
        }
        uniprotFile.createNewFile();
        FileWriter fstream = new FileWriter(uniprotFile);
        BufferedWriter out = new BufferedWriter(fstream);

        try (Transaction tx = graphDb.beginTx()) {
            ResourceIterator<Node> nodes = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            while(nodes.hasNext()){
                Node node = nodes.next();
                out.write(node.getProperty(PropertyType.DB_ID.toString()).toString() + "\n");
            }
            out.close();
            tx.success();
        }
        graphDb.shutdown();
    }

    /**
     * prints all phosphorylations in database except ptms from deleted UniProt IDs
     * @throws IOException
     */
    public void printAllPhosphorylations() throws IOException {
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        File uniprotFile = new File(outputFile + "/Phosphorylations.tsv");
        uniprotFile.createNewFile();

        FileWriter fstream = new FileWriter(uniprotFile);
        BufferedWriter out = new BufferedWriter(fstream);

        HashSet<String> ptms = new HashSet<>();

        try (Transaction tx = graphDb.beginTx()) {
            ResourceIterator<Node> phosNodes = graphDb.findNodes(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
            while (phosNodes.hasNext()) {
                Node phosNode = phosNodes.next();
                if(!phosNode.getProperty(PropertyType.LOCATION.toString()).equals("unknown")){
                    Iterable<Relationship> proteinRels = phosNode.getRelationships();
                    for(Relationship proteinRel : proteinRels){
                        Node protein = proteinRel.getEndNode();
                        Iterable<Relationship> UIDrelationships = protein.getRelationships(RelTypes.ID_BELONGS_TO);
                        for (Relationship UIDrelationship: UIDrelationships){
                            Node UID = UIDrelationship.getStartNode();
                            if(!UID.getProperty(PropertyType.STATUS.toString()).equals("Deleted")){
                                ptms.add(UID.getProperty(PropertyType.DB_ID.toString())
                                        +"_"+ phosNode.getProperty(PropertyType.TYPE.toString())
                                        + phosNode.getProperty(PropertyType.LOCATION.toString()));
                            }
                        }
                    }
                }

            }

            for(String ptm: ptms){
                out.write(ptm + "\n");
            }
            out.close();
            tx.success();
        }
        graphDb.shutdown();
    }

    /**
     * prints the number of things in the database with the given label
     * @param label
     */
    public void printLabelNumber(String label){
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> nodes = graphDb.findNodes(Label.label(label));
            Integer count = 0;
            for (ResourceIterator<Node> it = nodes; it.hasNext(); ) {
                Node node = it.next();
                count ++;
            }
            if(count.equals(0)){
                ResourceIterable<Label> allLabels = graphDb.getAllLabels();
                HashSet labels = new HashSet();
                for (Label l: allLabels) {
                    labels.add(l);
                }
                System.out.println("There are " +count+" "+label+"'s in this database: " + databaseDir);
                System.out.println("These are the labels available: ");
                System.out.println(labels);
            }else{
                System.out.println("There are " +count+" "+label+"'s in this database: " + databaseDir);
            }
        }

    }

    /**
     *  Prints a SIF file of the given database object
     *  Will also print an attribute file, linked by the node ID
     * @throws IOException
     */
    public void writeSIF() throws IOException {
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
        File SIF_FILE = new File(outputFile, "/SIF.sif");
        //File ATTR = new File(outputFile, "/attribute.tsv");

        FileWriter fstream = new FileWriter(outputFile + "/attribute.tsv");
        BufferedWriter ATTR = new BufferedWriter(fstream);

        try (Transaction tx = graphDb.beginTx()) {
            // Hashmap that keeps the printing order
            HashMap<Integer, String> orderKeeper = new HashMap<>(); //{1:DB_ID, 2: DisplayName}
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            Integer count = 0;
            for (String property: allPropertyKeys) {
                count++;
                orderKeeper.put(count, property);
            }

            String atrs = "Node_ID\t";
            List<Integer> properties2print = new ArrayList<>();
            for (Integer orderInt: orderKeeper.keySet()) {
                for (int j = 1; j < orderKeeper.keySet().size()+1; j++) {
                    if(orderInt == j){// in order
                        // get innerMap key
                        String property = orderKeeper.get(orderInt);
                            if(property.equals("DB_ID")){
                                atrs = atrs.concat( "Database_ID\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("DISPLAY_NAME")){
                                atrs = atrs.concat( "Display_Name\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("UNIPROT_ID")){
                                atrs = atrs.concat( "UniProt_ID\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("UNIPROT_NAME")){
                                atrs = atrs.concat( "UniProt_Gene_Name\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("TYPE")){
                                atrs = atrs.concat( "Type\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("LOCATION")){
                                atrs = atrs.concat( "Location\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("INTEGRATED")){
                                atrs = atrs.concat( "Integrated\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("STATUS")){
                                atrs = atrs.concat( "Status\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("KINASE")){
                                atrs = atrs.concat( "Kinase\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("TRANSCRIPTION_FACTOR")){
                                atrs = atrs.concat( "Transcription_Factor\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("CELL_SURFACE_RECEPTOR")){
                                atrs = atrs.concat( "Cell_Surface_Receptor\t");
                                properties2print.add(orderInt);
                            }else if(property.equals("DB_CONNECTION")){
                                atrs = atrs.concat( "Database_Link\t");
                                properties2print.add(orderInt);
                            }else if(property.contains("SUPPORT_SCORE_")){
                                atrs = atrs.concat( property+"\t");
                                properties2print.add(orderInt);
                            }else if(property.contains("ABUNDANCE_SCORE")){
                                atrs = atrs.concat( property+"\t");
                                properties2print.add(orderInt);
                            }
                    }
                }
            }

            atrs = atrs.concat("\n");
            ATTR.write(atrs);

            ResourceIterator<Node> peNodes = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            while (peNodes.hasNext()){
                Node pe = peNodes.next();

                atrs = pe.getId() + "\t";
                Map<String, Object> allProperties = pe.getAllProperties();
                for (Integer i: properties2print) {
                    String property = orderKeeper.get(i);
                    if(pe.hasProperty(property)){
                        atrs = atrs.concat(allProperties.get(property).toString() + "\t");
                    }else{
                        atrs = atrs.concat("\t");
                    }
                }
                atrs = atrs.concat("\n");
                ATTR.write(atrs);
            }
            ResourceIterator<Node> intNodes = graphDb.findNodes(Label.label(LabelTypes.INTERACTION.toString()));
            while (intNodes.hasNext()){
                Node pe = intNodes.next();

                atrs = pe.getId() + "\t";
                Map<String, Object> allProperties = pe.getAllProperties();
                for (Integer i: properties2print) {
                    String property = orderKeeper.get(i);
                    if(pe.hasProperty(property)){
                        atrs = atrs.concat(allProperties.get(property).toString() + "\t");
                    }else{
                        atrs = atrs.concat("\t");
                    }
                }
                atrs = atrs.concat("\n");
                ATTR.write(atrs);
            }
            ResourceIterator<Node> uidNodes = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            while (uidNodes.hasNext()){
                Node pe = uidNodes.next();

                atrs = pe.getId() + "\t";
                Map<String, Object> allProperties = pe.getAllProperties();
                for (Integer i: properties2print) {
                    String property = orderKeeper.get(i);
                    if(pe.hasProperty(property)){
                        atrs = atrs.concat(allProperties.get(property).toString() + "\t");
                    }else{
                        atrs = atrs.concat("\t");
                    }
                }
                atrs = atrs.concat("\n");
                ATTR.write(atrs);
            }
            ResourceIterator<Node> phosNodes = graphDb.findNodes(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
            while (phosNodes.hasNext()){
                Node pe = phosNodes.next();

                atrs = pe.getId() + "\t";
                Map<String, Object> allProperties = pe.getAllProperties();
                for (Integer i: properties2print) {
                    String property = orderKeeper.get(i);
                    if(pe.hasProperty(property)){
                        atrs = atrs.concat(allProperties.get(property).toString() + "\t");
                    }else{
                        atrs = atrs.concat("\t");
                    }
                }
                atrs = atrs.concat("\n");
                ATTR.write(atrs);
            }
            ResourceIterator<Node> modNodes = graphDb.findNodes(Label.label(LabelTypes.MODIFICATION.toString()));
            while (modNodes.hasNext()){
                Node pe = modNodes.next();

                atrs = pe.getId() + "\t";
                Map<String, Object> allProperties = pe.getAllProperties();
                for (Integer i: properties2print) {
                    String property = orderKeeper.get(i);
                    if(pe.hasProperty(property)){
                        atrs = atrs.concat(allProperties.get(property).toString() + "\t");
                    }else{
                        atrs = atrs.concat("\t");
                    }
                }
                atrs = atrs.concat("\n");
                ATTR.write(atrs);
            }
            ResourceIterator<Node> pathNodes = graphDb.findNodes(Label.label(LabelTypes.PATHWAY.toString()));
            while (pathNodes.hasNext()){
                Node pe = pathNodes.next();

                atrs = pe.getId() + "\t";
                Map<String, Object> allProperties = pe.getAllProperties();
                for (Integer i: properties2print) {
                    String property = orderKeeper.get(i);
                    if(pe.hasProperty(property)){
                        atrs = atrs.concat(allProperties.get(property).toString() + "\t");
                    }else{
                        atrs = atrs.concat("\t");
                    }
                }
                atrs = atrs.concat("\n");
                ATTR.write(atrs);
            }

            // SIF
            List<String> SIF = new ArrayList<String>();
            //SIF via all relationships
            ResourceIterable<Relationship> allRelationships = graphDb.getAllRelationships();
            for (Relationship relationship : allRelationships) {

                SIF.add(relationship.getStartNode().getId() +
                        "\t" + relationship.getType().toString() +
                        "\t" + relationship.getEndNode().getId() +
                        "\n");
            }

            SIF_FILE.createNewFile();
            // write SIF
            FileWriter writer = new FileWriter(SIF_FILE);
            for (String str : SIF) {
                writer.write(str);
            }

            tx.success();
            ATTR.close();
            writer.close();
        }
        graphDb.shutdown();
    }

    /**
     * Maps data onto database, giving proteins support score as well as abundance score
     * @param path2phosPeps
     * @throws IOException
     */
    public void mapMQPhosphopeps( File path2phosPeps, String mappingOption ) throws IOException{

        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        // read in uids from data
            // Ask for Leading Razor Protein Column
            // Ask for modified Peptide Column
            // Ask for column with abundance Score
            // Ask for name of Mapping

        Scanner scanner = new Scanner(System.in);
        System.out.println("\nPlease type the name of the column containing a single UniProt ID per row \n**There can be no spaces in the name**");
        String LRP_str = scanner.next();

        System.out.println("\nPlease type the name of the column containing the modified protein sequence\n**There can be no spaces in the name**");
        System.out.println("Example sequence formats: _(ac)AAAITDM(ox)ADLEELSRLS(ph)PLPPGS(ph)PGSAAR_  or AAAITDMADLEELSRLpSPLPPGpSPGSAAR ");
        String mod_pep = scanner.next();

        System.out.println("\nPlease type the name of the column containing the value you would like mapped\n**There can be no spaces in the name**");
        System.out.println("(e.g., the column that holds: p-value, SILAC_ratio, log2_intensity");
        String aVal = scanner.next();


        System.out.println("\nPlease type the name of the column you would like the mapped values associated with \n**There can be no spaces in the name**");
        System.out.println("(e.g., the column that holds:  stimulated, unstimulated, time_point_1");
        String experiment = scanner.next();

        /*
        String LRP_str = "Protein_ID";
        String mod_pep = "modified_peptide";
        String aVal = "log2Int";
        String experiment = "Time";

Protein_ID modified_peptide log2Int Time
         */

        /*
        String LRP_str = "UniProtID";
        String mod_pep = "mod_Seq";
        String aVal = "pVal";
        String experiment = "expr";
         */
        System.out.println("\nThank you, now extracting data from file");

        // read in human or mouse fasta file
        String species = getSpecies(graphDb);
        HashMap<String, String> fastaMap = new HashMap<>();
        if(species.equalsIgnoreCase("human")){
            // read in human fasta
            fastaMap = fasta2dict(humanUniProt);
        }else if(species.equalsIgnoreCase("mouse")){
            // read in mouse fasta
            fastaMap = fasta2dict(mouseUniProt);
        }

        //Create File to write inputs/Outputs to
        FileWriter fstream =  new FileWriter(outputFile + "/SupportingPeptides.tsv");
        BufferedWriter out = new BufferedWriter(fstream);
        out.write("UniProtID\tSupportingPeptideLineNumber:modification\n" );

        //Create File to write R inputs to
        FileWriter fstream1 =  new FileWriter(outputFile + "/R_Mapping_input.tsv");
        BufferedWriter out1 = new BufferedWriter(fstream1);
        out1.write("Experiment\tpepMap\tpepNotMap\tUIDMap\tUIDnotMap\tphosdTo\tnotPhosdTo\tpfTo\tnotPFTo\tcplxTo\tnotCplxTo\tintTo\tnotIntTo\n");


        // Read in data into Dict
        BufferedReader BR = new BufferedReader(new FileReader(path2phosPeps));
        String line = "";
        int ID = 0;


        Integer leadingProtCol =0;
        Integer modPepCol = 0;
        Integer aScoreCol = 0;
        Integer experimentCol = 0;

        Boolean lpBool= true;
        Boolean modPepBool= true;
        Boolean aSBool= true;
        Boolean expBool= true;

        // first get all experiemnt names by reading through the file
        HashSet<String> experiments = new HashSet<>();
        while ((line = BR.readLine()) != null) {
            ID ++;
            String[] columns = line.split("\t");
            List<String> cols = new ArrayList<>();

            Integer colLen = columns.length;
            for (int i = 0; i < colLen; i++) {
                cols.add(columns[i]);
                if(columns[i].equals(experiment)){
                    experimentCol = i;
                    expBool = false;
                }
                else if (columns[i].equalsIgnoreCase(LRP_str)){
                    leadingProtCol = i;
                    lpBool = false;
                }
                else if(columns[i].equals(mod_pep)){
                    modPepCol = i;
                    modPepBool = false;
                }
                else if(columns[i].equals(aVal)){
                    aScoreCol = i;
                    aSBool = false;
                }
            }

            // Make sure the columns input exist
            if(lpBool){
                throw new IOException("There is no "+LRP_str+" column in the input data file" +
                        "\nThese are the columns available \n" + cols);

            }else if(modPepBool){
                throw new IOException("There is no "+mod_pep+" column in the input data file" +
                        "\nThese are the columns available \n" + cols);
            }else if(aSBool){
                throw new IOException("There is no "+aVal+" column in the input data file" +
                        "\nThese are the columns available \n" + cols);
            }else if(expBool){
                throw new IOException("There is no "+experiment+" column in the input data file" +
                        "\nThese are the columns available \n" + cols);
            }

            if (!columns[experimentCol].equals(experiment)){
                experiments.add(columns[experimentCol]);
            }
        }
        System.out.println("Experiments to map"+ experiments);

        HashMap<Node, HashSet<Node>> pathSets = getPathSets(graphDb);

        // next for each different experiment name, only gather the data in those experiments.
        for (String exprName:experiments) {
            //Create File to write inputs/Outputs to
            FileWriter fstream2 =  new FileWriter(outputFile + "/PhosphoMappingReport_"+exprName+".tsv");
            BufferedWriter out2 = new BufferedWriter(fstream2);

            System.out.println("\nCurrently Mapping: "+ exprName);
            String supportScoreString = "SUPPORT_SCORE_" + exprName;
            String abundanceScoreString = "ABUNDANCE_SCORE_" + exprName;
            String scoredByString = "SCORED_BY" + exprName;
            String mappedString = "MAPPED_" + exprName;

            Double maxAbundance = Double.NEGATIVE_INFINITY;
            Double minAbundance = Double.POSITIVE_INFINITY;


            // extract:
            HashMap<String, HashSet<String>> ID2LRP = new HashMap<>(); // {LRP: id1, id2, id3}
            HashMap<String, Double> ID2eScore = new HashMap<>();  // {ID: p-value}
            // generate from input:
            HashMap<String, HashSet<String>> ID2TopMods = new HashMap<>(); // {ID: S_19,T_18, Y_42}
            HashMap<String, String> ID2Start = new HashMap<>(); // {ID: start of pep}
            HashMap<String, String> ID2End = new HashMap<>(); // {ID: end of pep}

            BR = new BufferedReader(new FileReader(path2phosPeps));
            line = "";
            ID = 0;


            // get columns of interest from data file
            while ((line = BR.readLine()) != null) {
                ID ++;
                String[] columns = line.split("\t");

                if(columns[experimentCol].equals(exprName)){ // only per experiment
                    // extract LRP's
                    if(ID2LRP.containsKey(columns[leadingProtCol])){
                        HashSet<String> temp = ID2LRP.get(columns[leadingProtCol]);
                        temp.add(String.valueOf(ID));
                        ID2LRP.put(columns[leadingProtCol], temp);
                    }else{
                        HashSet<String> temp = new HashSet<>();
                        temp.add(String.valueOf(ID));
                        ID2LRP.put(columns[leadingProtCol], temp);
                    }

                    //extract aScores:
                    if(columns[aScoreCol].isEmpty() | columns[aScoreCol].equalsIgnoreCase(aVal)) {
                        ID2eScore.put(String.valueOf(ID), null);
                    }else if(columns[aScoreCol].equalsIgnoreCase("NA")){
                        ID2eScore.put(String.valueOf(ID), null);
                    }else{
                        ID2eScore.put(String.valueOf(ID), Double.valueOf(columns[aScoreCol]));
                    }

                    if(!columns[modPepCol].equalsIgnoreCase(mod_pep)){
                        // generate Start
                        String start = generateStartOfPep(fastaMap, columns[leadingProtCol],columns[modPepCol]);
                        if(start.isEmpty()){
                            continue;
                        }

                        ID2Start.put(String.valueOf(ID), start);


                        // generate End
                        String end = generateEndOfPep(start, columns[modPepCol]);
                        ID2End.put(String.valueOf(ID), end);

                        // generate Top Mods
                        HashSet<String> mods = generateModList(ID,columns[modPepCol], ID2Start);
                        ID2TopMods.put(String.valueOf(ID), mods);
                    }
                }
            }

/*
            System.out.println("LRP: "+ ID2LRP);
            System.out.println("ascore"+ ID2eScore);
            System.out.println("id2start: " + ID2Start);
            System.out.println("id2end: " + ID2End);
            System.out.println("modlist: " + ID2TopMods);
*/


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

                        UID.setProperty(mappedString, "Mapped");

                        // GET UID
                        Iterable<Relationship> protRelationships = UID.getRelationships(RelTypes.ID_BELONGS_TO);

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


                        for (Relationship protRltnshp : protRelationships) {
                            //System.out.println("//////////////////////////////////////////////////////////////////////////////");
                            double proteoformScore = 0.0;
                            Node protein = protRltnshp.getEndNode();

                            // GET PF MODs
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

                            //////////////////////////////////////////////////////////////////////////////////////////
                            // setting the support score
                            // for each peptide
                            HashMap<String, Double> pepID2sScore = new HashMap<>();
                            for(String PeptideID: PeptideIDs) {
                                if(!ID2Start.containsKey(PeptideID)){
                                    continue;
                                }
                                HashSet<String> pepMods = ID2TopMods.get(PeptideID);
                                Integer start = Integer.valueOf(ID2Start.get(PeptideID));
                                Integer end = Integer.valueOf(ID2End.get(PeptideID));

                                out.write(LRP + "\t" + PeptideID + ":" + pepMods +"\n");

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

                                /*
                                System.out.println();
                                System.out.println();
                                //System.out.println(LRP);
                                System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()));
                                Integer pepid = Integer.valueOf(PeptideID);
                                System.out.println(pepid);
                                System.out.println("Peptide Mods "+ pepid+": " + pepMods);
                                System.out.println("POI in peptide range("+start+" , "+end+") : "+ POIwithinPep);
                                System.out.println("All mods in pf: " + proteoformModList);
                                 */

                                double unModPeptideScore = 0.0;
                                double peptideScore = 0.0;
                                if(POIwithinPep.isEmpty()){ // no mods on uid
                                    if(!pepMods.isEmpty()){ // if mods in pep
                                        //System.out.println("Pep covers non modified region but has mods itself += 0.9");
                                        peptideScore+= 0.9;
                                    }else{ // no mods in pep (perf score without multiplying by 1.5)
                                        //System.out.println("Pep covers non modified region and has no mods += 1");
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

                                    // remove all POI's from copy list so your left with mods not seen in the db
                                    pepModCopy.removeAll(UIDmodList);

                                    if(!pepModCopy.isEmpty()){ // if pep has mods that are not POI, not exact match changed * 0.9
                                        peptideScore = peptideScore*0.9;
                                        //System.out.println("Modified pep had additional mods (not in db) *.0.9");
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
                            //System.out.println("pf score: " + proteoformScore);


                            //System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()) +" "+pepID2sScore);

                            // If protein already has a score - must have multi UIDs
                            // so take the uid that gives it the highest score
                            //System.out.println( "PF score: "+ proteoformScore);
                            if(protein.hasProperty(supportScoreString)){
                                double currentScore = Double.parseDouble(protein.getProperty(supportScoreString).toString());
                                if(proteoformScore > currentScore){
                                    protein.setProperty(supportScoreString, (Math.round(proteoformScore*100.00))/100.00);
                                    protein.setProperty(scoredByString, LRP);
                                }
                            }else{ // no previous score - set it
                                protein.setProperty(supportScoreString, (Math.round(proteoformScore*100.00))/100.00);
                            }

                            //System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()) +" pf sScore:"+ protein.getProperty(supportScoreString));

                            pepsSpanningModTemp = pepsSpanningModTemp/getLength(protRelationships);
                            pepsSpanningMod += pepsSpanningModTemp;
                            pepsSpanningModTemp =0;

                            //////////////////////////////////////////////////////////////////////////////////////////
                            // Setting abundance score:
                            // find the Peptide ID's with the biggest score

                            HashMap<String, Double> maxMin = new HashMap<>();
                            if(mappingOption.equalsIgnoreCase("HighestSupport")){
                                maxMin = abundanceScoreAvgHighestSupport(graphDb, pepID2sScore, ID2eScore, protein, abundanceScoreString, maxAbundance, minAbundance);
                            }else if(mappingOption.equalsIgnoreCase("Max")){
                                maxMin = abundanceScoreHighestOfAllPeps(graphDb, PeptideIDs, ID2eScore, protein, abundanceScoreString, maxAbundance, minAbundance);
                            }else if (mappingOption.equalsIgnoreCase("Mean")){
                                maxMin = abundanceScoreMean(graphDb, PeptideIDs, ID2eScore, protein, abundanceScoreString, maxAbundance, minAbundance);
                            }else if(mappingOption.equalsIgnoreCase("Median")){
                                maxMin = abundanceScoreMedian(graphDb, PeptideIDs, ID2eScore, protein, abundanceScoreString, maxAbundance, minAbundance);
                            }else if(mappingOption.equalsIgnoreCase("Extreme")){
                                maxMin = abundanceScoreMostExtremeOfAllPeps(graphDb, PeptideIDs, ID2eScore, protein, abundanceScoreString, maxAbundance, minAbundance);
                            }else{
                                throw new InputException("Abundance Score Mapping options are \"HighestSupport\", \"Max\",\"Median\",\"Mean\".");
                            }
                            //System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()) +" pf aScore:"+ protein.getProperty(abundanceScoreString));

                            maxAbundance = maxMin.get("Max");
                            minAbundance = maxMin.get("Min");

                            maxAbundance = (Math.round(maxAbundance*100.0))/100.0;
                            minAbundance = (Math.round(minAbundance*100.0))/100.0;

                            //////////////////////////////////////////////////////////////////////////////////////////

                        }
                    }
                }

                out2.write("Data Statistics:");
                out2.write("\nThe number of peptides mapped: " + pepsMatched);
                out2.write("\nThe number of peptides not mapped: " + (ID -pepsMatched) );
                out2.write("\nThe number of peptides spanning a modified region of a proteoform: " + Math.round(pepsSpanningMod));
                out2.write("\nThe number of UniProt ID's in data that mapped to database: " + UIDsMatched);
                out2.write("\nThe number of UniProt ID's in data that did not map to database: " + UIDsNotMatched); //TODO WRONG? - ONLY COUNTS UIDS IN DB
                out2.write("\nThe range of abundances: " + minAbundance +" to "+ maxAbundance);

                out1.write(exprName +"\t"+pepsMatched+"\t"+(ID -pepsMatched) +"\t"+UIDsMatched +"\t"+UIDsNotMatched +"\t");

                tx.success();
            }

            mapComplexs(graphDb, supportScoreString, abundanceScoreString,  scoredByString, mappedString);
            addRelWeights(graphDb, supportScoreString, maxAbundance, minAbundance);
            mappingReport(out2, graphDb, supportScoreString, out1, pathSets);
            out2.close();

        }// expr for loop
        System.out.println("\n\nPlease enter the following line into your console to generate diagnostic plots");
        System.out.println("Rscript -e \"rmarkdown::render('proportionPlots.Rmd')\"");
        graphDb.shutdown();
        out.close();
        out1.close();
    }

    HashMap<String, Double> abundanceScoreAvgHighestSupport(GraphDatabaseService graphDb,
                                                            HashMap<String, Double> pepID2sScore, // {pepID: supportScore}
                                                            HashMap<String, Double> ID2eScore, //{pepID: abundanceScore}
                                                            Node protein, //current protein score is being set for
                                                            String abundanceScoreString, // the cirrent experiment the score is being set for
                                                            Double maxAbundance, // max abundance score in whole ntwk
                                                            Double minAbundance){  // min abundance score in whole ntwk
        HashMap<String, Double> abundanceScores = new HashMap<>();

        try(Transaction tx = graphDb.beginTx()){

            // get the pepId with the largest support score
            Map.Entry<String, Double> maxEntry = null;
            for (Map.Entry<String, Double> entry: pepID2sScore.entrySet()){
                if(maxEntry == null || entry.getValue().compareTo(maxEntry.getValue()) > 0 ){ // returns >0 if the value we are looking at is > than the current maxEntry
                    maxEntry = entry;
                }
            }

            // may be multi pep's with largest score, add all to hashset
            HashSet<String> maxPeps = new HashSet<>();
            for(String pep: pepID2sScore.keySet()){
                Double score = pepID2sScore.get(pep); // if multi with max score
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
                    if(ID2eScore.get(pep) != null){
                        eScores[pepcount] = ID2eScore.get(pep);
                    }
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
                    eScore = ID2eScore.get(pep);
                }
            }

            // If protein already has a score - must have multi UIDs
            // so take the uid that gives it the highest score
            //System.out.println( "PF score: "+ proteoformScore);
            if(protein.hasProperty(abundanceScoreString)){
                double currentScore = Double.parseDouble(protein.getProperty(abundanceScoreString).toString());
                if(eScore > currentScore){
                    protein.setProperty(abundanceScoreString, (Math.round(eScore*100.00))/100.00);
                }
            }else{ // no previous score - set it
                protein.setProperty(abundanceScoreString, (Math.round(eScore*100.00))/100.00);
            }

            //System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()) +" pf aScore:"+ protein.getProperty(abundanceScoreString));

            // get max and min Scores
            if(eScore > maxAbundance){
                maxAbundance = eScore;
            }

            if(eScore < minAbundance){
                minAbundance = eScore;
            }
            abundanceScores.put("Max", maxAbundance);
            abundanceScores.put("Min", minAbundance);

            tx.success();
        }

        return abundanceScores;
    }

    HashMap<String, Double> abundanceScoreHighestOfAllPeps(GraphDatabaseService graphDb,
                                                                   HashSet<String> PeptideIDs, // set of Ids mapping to this proteofrom
                                                                   HashMap<String, Double> ID2eScore, //{pepID: abundanceScore}
                                                                   Node protein, //current protein score is being set for
                                                                   String abundanceScoreString, // the cirrent experiment the score is being set for
                                                                   Double maxAbundance, // max abundance score in whole ntwk
                                                                   Double minAbundance){  // min abundance score in whole ntwk

        HashMap<String, Double> abundanceScores = new HashMap<>();

        try(Transaction tx = graphDb.beginTx()){

            // get a map of id's to abundances for just this proteoform
            HashMap<String, Double> pepID2eScore = new HashMap<>();
            for(String pepID: PeptideIDs){
                pepID2eScore.put(pepID, ID2eScore.get(pepID));
            }

            // get the peptide with the max escore
            Map.Entry<String, Double> maxEntry = null;
            for (Map.Entry<String, Double> entry: pepID2eScore.entrySet()){
                if(maxEntry == null || entry.getValue().compareTo(maxEntry.getValue()) > 0 ){ // returns >0 if the value we are looking at is > than the current maxEntry
                    maxEntry = entry;
                }
            }
            // now you have the peptide Id with the largest abundance in maxEntry
            Double eScore = maxEntry.getValue();

            // If protein already has a score - must have multi UIDs
            // so take the uid that gives it the highest score
            //System.out.println( "PF score: "+ proteoformScore);
            if(protein.hasProperty(abundanceScoreString)){
                double currentScore = Double.parseDouble(protein.getProperty(abundanceScoreString).toString());
                if(eScore > currentScore){
                    protein.setProperty(abundanceScoreString, (Math.round(eScore*100.00))/100.00);
                }
            }else{ // no previous score - set it
                protein.setProperty(abundanceScoreString, (Math.round(eScore*100.00))/100.00);
            }


            // get max and min Scores
            if(eScore > maxAbundance){
                maxAbundance = eScore;
            }

            if(eScore < minAbundance){
                minAbundance = eScore;
            }

            abundanceScores.put("Max", maxAbundance);
            abundanceScores.put("Min", minAbundance);


            tx.success();
        }
        return abundanceScores;
    }
    
    HashMap<String, Double> abundanceScoreMostExtremeOfAllPeps(GraphDatabaseService graphDb,
                                                           HashSet<String> PeptideIDs, // set of Ids mapping to this proteofrom
                                                           HashMap<String, Double> ID2eScore, //{pepID: abundanceScore}
                                                           Node protein, //current protein score is being set for
                                                           String abundanceScoreString, // the cirrent experiment the score is being set for
                                                           Double maxAbundance, // max abundance score in whole ntwk
                                                           Double minAbundance){  // min abundance score in whole ntwk

        HashMap<String, Double> abundanceScores = new HashMap<>();

        try(Transaction tx = graphDb.beginTx()){

            // get a map of id's to abundances for just this proteoform
            HashMap<String, Double> pepID2eScore = new HashMap<>();
            for(String pepID: PeptideIDs){
                pepID2eScore.put(pepID, ID2eScore.get(pepID));
            }

            // get the peptide with the max escore
            Map.Entry<String, Double> maxEntry = null;
            for (Map.Entry<String, Double> entry: pepID2eScore.entrySet()){
                if(maxEntry == null || entry.getValue().compareTo(Math.abs(maxEntry.getValue())) > 0 ){ // returns >0 if the |value| we are looking at is > than the current maxEntry
                    maxEntry = entry;
                }
            }
            // now you have the peptide Id with the largest abundance in maxEntry
            Double eScore = maxEntry.getValue();

            // If protein already has a score - must have multi UIDs
            // so take the uid that gives it the highest score
            //System.out.println( "PF score: "+ proteoformScore);
            if(protein.hasProperty(abundanceScoreString)){
                double currentScore = Double.parseDouble(protein.getProperty(abundanceScoreString).toString());
                if(eScore > currentScore){
                    protein.setProperty(abundanceScoreString, (Math.round(eScore*100.00))/100.00);
                }
            }else{ // no previous score - set it
                protein.setProperty(abundanceScoreString, (Math.round(eScore*100.00))/100.00);
            }


            // get max and min Scores
            if(eScore > maxAbundance){
                maxAbundance = eScore;
            }

            if(eScore < minAbundance){
                minAbundance = eScore;
            }

            abundanceScores.put("Max", maxAbundance);
            abundanceScores.put("Min", minAbundance);


            tx.success();
        }
        return abundanceScores;
    }

    HashMap<String, Double>  abundanceScoreMean(GraphDatabaseService graphDb,
                                                        HashSet<String> PeptideIDs, // set of Ids mapping to this proteofrom
                                                        HashMap<String, Double> ID2eScore, //{pepID: abundanceScore}
                                                        Node protein, //current protein score is being set for
                                                        String abundanceScoreString, // the cirrent experiment the score is being set for
                                                        Double maxAbundance, // max abundance score in whole ntwk
                                                        Double minAbundance){  // min abundance score in whole ntwk

        HashMap<String, Double> abundanceScores = new HashMap<>();

        try(Transaction tx = graphDb.beginTx()){
            // get a map of id's to abundances for just this proteoform
            HashMap<String, Double> pepID2eScore = new HashMap<>();
            for(String pepID: PeptideIDs){
                pepID2eScore.put(pepID, ID2eScore.get(pepID));
            }

            //get the avg for all abundances
            Collection<Double> values = pepID2eScore.values();
            Double total = 0.0;
            for (Double i: values){
                total += i;
            }
            Double avg = total/values.size();

            // If protein already has a score - must have multi UIDs
            // so take the uid that gives it the highest score
            //System.out.println( "PF score: "+ proteoformScore);
            if(protein.hasProperty(abundanceScoreString)){
                double currentScore = Double.parseDouble(protein.getProperty(abundanceScoreString).toString());
                if(avg > currentScore){
                    protein.setProperty(abundanceScoreString, (Math.round(avg*1000d))/1000d);
                }
            }else{ // no previous score - set it
                protein.setProperty(abundanceScoreString, (Math.round(avg*1000d))/1000d);
            }

            // get max and min Scores
            if(avg > maxAbundance){
                maxAbundance = avg;
            }

            if(avg < minAbundance){
                minAbundance = avg;
            }

            abundanceScores.put("Max", maxAbundance);
            abundanceScores.put("Min", minAbundance);

            tx.success();
        }

        return abundanceScores;

    }

    HashMap<String, Double> abundanceScoreMedian(GraphDatabaseService graphDb,
                                                         HashSet<String> PeptideIDs, // set of Ids mapping to this proteofrom
                                                         HashMap<String, Double> ID2eScore, //{pepID: abundanceScore}
                                                         Node protein, //current protein score is being set for
                                                         String abundanceScoreString, // the cirrent experiment the score is being set for
                                                         Double maxAbundance, // max abundance score in whole ntwk
                                                         Double minAbundance){  // min abundance score in whole ntwk

        HashMap<String, Double> abundanceScores = new HashMap<>();

        try(Transaction tx = graphDb.beginTx()){

            HashMap<String, Double> pepID2eScore = new HashMap<>();
            for(String pepID: PeptideIDs){
                pepID2eScore.put(pepID, ID2eScore.get(pepID));
            }

            //get the median for all abundances
            Collection<Double> values = pepID2eScore.values();
            ArrayList<Double> vals = new ArrayList<>(values);
            Collections.sort(vals);

            Double median = 0.0;
            if(vals.size()%2 ==0){ // if even
                median = ((vals.get((vals.size()/ 2)-1) + vals.get(vals.size()/2))/2);
            }else{ // if odd
                median = vals.get(((vals.size()+1)/2)-1);
            }


            // If protein already has a score - must have multi UIDs
            // so take the uid that gives it the highest score
            //System.out.println( "PF score: "+ proteoformScore);
            if(protein.hasProperty(abundanceScoreString)){
                double currentScore = Double.parseDouble(protein.getProperty(abundanceScoreString).toString());
                if(median > currentScore){
                    protein.setProperty(abundanceScoreString, (Math.round(median*1000d))/1000d);
                }
            }else{ // no previous score - set it
                protein.setProperty(abundanceScoreString, (Math.round(median*1000d))/1000d);
            }

            // get max and min Scores
            if(median > maxAbundance){
                maxAbundance = median;
            }

            if(median < minAbundance){
                minAbundance = median;
            }

            abundanceScores.put("Max", maxAbundance);
            abundanceScores.put("Min", minAbundance);

            tx.success();
        }

        return abundanceScores;

    }

    /**
     * A function that reports statistics for mapping
     * @param out2
     * @param graphDb
     * @param supportScoreString
     * @throws IOException
     */
    private void mappingReport(BufferedWriter out2,
                               GraphDatabaseService graphDb,
                               String supportScoreString,
                               BufferedWriter out1,
                               HashMap<Node, HashSet<Node>> pathSet) throws IOException {
        try (Transaction tx = graphDb.beginTx()) {
            // Phosn stats
            Integer phosnsMappedTo = 0;
            Integer phosnsNotMappedTo = 0;
            ResourceIterator<Node> phosns = graphDb.findNodes(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
            while(phosns.hasNext()){
                Node phos = phosns.next();
                Iterable<Relationship> phosnrelationships = phos.getRelationships(RelTypes.PHOSPHORYLATION);
                for (Relationship phosnRel: phosnrelationships){
                    Node protein = phosnRel.getEndNode();
                    if(protein.hasProperty(supportScoreString)){
                        phosnsMappedTo++;
                    }else{
                        phosnsNotMappedTo++;
                    }
                }
            }

            // pf stats
            Integer pfsMappedTo = 0;
            Integer pfsNotMappedTo = 0;
            Integer kinsMappedTo = 0;
            Integer kinsNotMappedTo = 0;
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while (proteins.hasNext()){
                Node protein = proteins.next();
                if(protein.hasProperty(supportScoreString)){
                    pfsMappedTo++;
                }else{
                    pfsNotMappedTo++;
                }
                if (protein.hasProperty(PropertyType.KINASE.toString())){
                    if(protein.hasProperty(PropertyType.KINASE.toString())){
                        kinsMappedTo++;
                    }else{
                        kinsNotMappedTo++;
                    }
                }
            }

            // cplx stats
            Integer cplxsMappedTo = 0;
            Integer cplxsNotMappedTo = 0;
            ResourceIterator<Node> complexes = graphDb.findNodes(Label.label("Complex"));
            while(complexes.hasNext()){
                Node complex = complexes.next();
                if (complex.hasProperty(supportScoreString)){
                    cplxsMappedTo++;
                }else{
                    cplxsNotMappedTo++;
                }
            }

            // UID stats
            Integer uidsMappedTo = 0;
            Integer uidsNotMappedTo = 0;
            ResourceIterator<Node> uids = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            while (uids.hasNext()){
                Node uid = uids.next();
                Boolean measured = false;
                for (Relationship relationship : uid.getRelationships()) {
                    Node prot = relationship.getEndNode();
                    if (prot.hasProperty(supportScoreString)){
                        measured = true;
                    }
                }
                if (measured){
                    uidsMappedTo++;
                }else{
                    uidsNotMappedTo++;
                }
            }

            // Integrated Stats
            Integer intsMappedTo = 0;
            Integer intsNotMappedTo = 0;
            ResourceIterator<Node> integrateds = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()), PropertyType.INTEGRATED.toString(), "true");
            while (integrateds.hasNext()){
                Node integrated = integrateds.next();
                if (integrated.hasLabel(Label.label("Protein"))){
                    if (integrated.hasProperty(supportScoreString)){
                        intsMappedTo++;
                    }else{
                        intsNotMappedTo++;
                    }
                }
            }

            out2.write("\n\nNetwork statistics: ");
            out2.write("\nThe number of phosporylations in the database that were mapped to: " + phosnsMappedTo);
            out2.write("\nThe number of phosporylations in the database that were not mapped to: " + phosnsNotMappedTo);
            out2.write("\nThe number of proteoforms in the database that were mapped to: " + pfsMappedTo);
            out2.write("\nThe number of proteoforms in the database that were not mapped to: " + pfsNotMappedTo);
            out2.write("\nThe number of kinases in the database that were mapped to: " + kinsMappedTo);
            out2.write("\nThe number of kinases in the database that were not mapped to: " + kinsNotMappedTo);
            out2.write("\nThe number of complexes in the database that were mapped to: " + cplxsMappedTo);
            out2.write("\nThe number of complexes in the database that were not mapped to: " + cplxsNotMappedTo);
            out2.write("\nThe number of UniProt IDs in the database that were mapped to: " + uidsMappedTo);
            out2.write("\nThe number of UniProt IDs in the database that were not mapped to: " + uidsNotMappedTo);
            out2.write("\nThe number of integrated nodes from PhosphoSitePlus in the database that were mapped to: " + intsMappedTo);
            out2.write("\nThe number of integrated nodes from PhosphoSitePlus in the database that were not mapped to: " + intsNotMappedTo);

            out1.write(phosnsMappedTo+"\t"+phosnsNotMappedTo+"\t"+pfsMappedTo+"\t"+pfsNotMappedTo+"\t"+cplxsMappedTo+"\t"+cplxsNotMappedTo+"\t"+intsMappedTo+"\t"+intsNotMappedTo+"\n");



            // Cell location Stats
            ResourceIterator<Node> nodes = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            HashMap<String, Integer> cellLocNums = new HashMap<>();
            while (nodes.hasNext()){
                Node node = nodes.next();
                if (node.hasProperty(PropertyType.LOCATION.toString()) & node.hasProperty(supportScoreString)){
                    String location = node.getProperty(PropertyType.LOCATION.toString()).toString();
                    if(cellLocNums.containsKey(location)){
                        Integer integer = cellLocNums.get(location);
                        integer++;
                        cellLocNums.put(location, integer);
                    }else{
                        cellLocNums.put(location, 1);
                    }
                }
            }

            // create a set of the numbers of things measured (to be deleted from)
            HashSet<Integer> uniqueNumMeasuredSet = new HashSet<>(cellLocNums.values());

            // now you have 1 map CellLoc:#measured
            // print largest, delete largest from list
            Collection<Integer> values = cellLocNums.values();
            List<Integer> numMeasured = new ArrayList<>(values);

            out2.write("\n\nNumber of proteins and complexes mapped in each cellular location:");
            if(uniqueNumMeasuredSet.size() > 1){
                Integer max = Collections.max(numMeasured);
                for (int i = 0; i < uniqueNumMeasuredSet.size(); i++) {
                    for (String key: cellLocNums.keySet()) {
                        if (cellLocNums.get(key).equals(max)){
                            out2.write("\n"+key+": "+ max);
                        }
                    }
                    numMeasured.remove(max);
                    max = Collections.max(numMeasured);
                }
            }

            tx.success();
        }

        // next get # things measured in path Sets
        HashMap<Node, Double> proportionMeasured = getPathProportions(graphDb,pathSet, supportScoreString);

        // print in order?
        Collection<Double> values1 = proportionMeasured.values();
        List<Double> proportions = new ArrayList<>(values1);
        Collections.sort(proportions);
        Collections.reverse(proportions);
        //List<Post> myLastPosts = posts.subList(posts.size()-40, posts.size());
        out2.write("\n\nTop Pathways with the highest proportions of molecules mapped: ");
        out2.write("\nPathways are only displayed if they contain > 10 mappable molecules: ");
        if(values1.size() > 100){
            Double largestProp = (Double) Collections.max(proportions);
            for (int i = 0; i < 100; i++) {
                for(Node key: proportionMeasured.keySet()){
                    if(proportionMeasured.get(key) == largestProp & pathSet.get(key).size() > 10){
                        try(Transaction tx = graphDb.beginTx()){
                            out2.write("\n");
                            recursePrintSubPaths( key, proportionMeasured, pathSet, -1, out2);
                            tx.success();
                        }
                    }
                }
                proportions.remove(largestProp);
                largestProp = (Double) Collections.max(proportions);
            }
        }

    }

    private void recursePrintSubPaths(  Node pathNode, HashMap<Node, Double> propMeasured,HashMap<Node, HashSet<Node>> pathSet, Integer count, BufferedWriter out) throws IOException {

            Iterable<Relationship> subPathRels = pathNode.getRelationships(RelTypes.SUB_PATHWAY, Direction.OUTGOING);
            if(getLength(subPathRels).equals(0)){
                count ++;
                out.write("\n");
                for (int i = 0; i < count; i++) {
                    out.write("\t");
                }
                out.write(pathNode.getProperty(PropertyType.DISPLAY_NAME.toString()) +": "+ propMeasured.get(pathNode) + " (Size: " + pathSet.get(pathNode).size()+ ")");
                return;
            }else{
                out.write("\n");
                count ++;
                for (int i = 0; i < count; i++) {
                    out.write("\t");
                }
                out.write(pathNode.getProperty(PropertyType.DISPLAY_NAME.toString())+": "+ propMeasured.get(pathNode)+ " (Size: " + pathSet.get(pathNode).size()+ ")");

                for (Relationship subPathRel: subPathRels) {
                    recursePrintSubPaths( subPathRel.getEndNode(),propMeasured, pathSet,  count, out);
                }
            }

    }

    private HashMap<Node, HashSet<Node>> getPathSets(GraphDatabaseService graphDb){
        HashMap<Node, HashSet<Node>> pathMap = new HashMap<>();

        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> allPathNodes = graphDb.findNodes(Label.label(LabelTypes.PATHWAY.toString()));
            // find all root nodes
            HashSet<Node> roots = new HashSet<>();
            while(allPathNodes.hasNext()){
                Node pathNode = allPathNodes.next();
                pathMap.put(pathNode, new HashSet<>());
                if (!pathNode.hasRelationship(RelTypes.SUB_PATHWAY, Direction.INCOMING)){
                    roots.add(pathNode);
                }
            }

            // keep track of depth via reverse BFS
            TraversalDescription rbfs = graphDb.traversalDescription()
                    .order(BranchOrderingPolicies.POSTORDER_BREADTH_FIRST)
                    .relationships(RelTypes.SUB_PATHWAY, Direction.OUTGOING)
                    .uniqueness(Uniqueness.NONE);

            // for each root node bfs up, incrementally adding to sets
            for (Node pathNode: roots) {
                Traverser traverse = rbfs.traverse(pathNode);
                for(Path pathPath: traverse){
                    HashSet<Node> leaves = new HashSet<>();
                    // for the final node in each path, iterate all of it's components
                    for (Relationship relationship : pathPath.endNode().getRelationships(RelTypes.PATHWAY_COMPONENT, Direction.OUTGOING)) {
                        //if the component is a PE, add to set
                        if(relationship.getEndNode().hasLabel(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()))){
                            leaves.add(relationship.getEndNode());
                        }
                    }
                    HashSet<Node> nodes = pathMap.get(pathPath.endNode());
                    nodes.addAll(leaves);
                    pathMap.put(pathPath.endNode(), nodes);

                    if(pathPath.length() > 0){
                        HashSet<Node> nodesParent = pathMap.get(pathPath.lastRelationship().getStartNode());
                        nodesParent.addAll(nodes);
                        pathMap.put(pathPath.lastRelationship().getStartNode(), nodesParent);
                    }
                    //System.out.println(pathPath);
                }
            }
            tx.success();
        }
        return pathMap;
    }

    private HashMap<Node, Double> getPathProportions(GraphDatabaseService graphDb, HashMap<Node, HashSet<Node>> pathSets, String supportScoreStr){
        HashMap<Node, Double> pathProp = new HashMap<>();

        try(Transaction tx = graphDb.beginTx()){
            for(Map.Entry<Node, HashSet<Node>> entry: pathSets.entrySet()){
                HashSet<Node> pathNodes = entry.getValue();
                Double numMeasured = 0.0;
                for (Node n: pathNodes) {
                    if(n.hasProperty(supportScoreStr)){
                        numMeasured ++;
                    }
                }
                Double prop = numMeasured/pathNodes.size();
                //System.out.println(entry.getKey() + ": "+ numMeasured +"+"+pathNodes.size() + "= "+ prop + " -> " + pathSets.get(entry.getKey()));
                pathProp.put(entry.getKey(), prop);
            }

            tx.success();
        }
        return pathProp;
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

    /**
     * Takes in a fasta file and returns a hashmap of ID:sequence
     * @param url
     * @return
     */
    HashMap<String, String> fasta2dict(String url)  {

        HashMap<String, String> fasta = new HashMap<>();
        URL fasta_url = null;
        try {
            fasta_url = new URL(url);
        } catch (MalformedURLException e) {
            e.printStackTrace();
            System.exit(1);
        }

        try  {
            BufferedReader BR = new BufferedReader(new InputStreamReader(fasta_url.openStream()));

            boolean first = true;
            String line;

            String UID = "";
            while((line = BR.readLine()) != null){
                line = line.trim();

                if(line.charAt(0) == '>'){
                    if (first){
                        first = false;
                    }
                    // need to extract UID
                    Pattern p = Pattern.compile(UID_PATTERN);
                    Matcher m = p.matcher(line);
                    m.find();
                    UID = m.group(0);
                    fasta.put(UID, "");
                }else{
                    String seqFrag = fasta.get(UID);
                    seqFrag = seqFrag.concat(line);
                    fasta.put(UID, seqFrag);
                }
            }
        }catch (IOException e){
            e.printStackTrace();
            System.err.println("Could not load file: " + url);
            System.exit(1);
        }

        return fasta;

    }

    /**
     * adds support scores and enrchment scores to Complexes
     * @param graphDb
     * @param supportScoreString
     * @param abundanceScoreString
     * @param scoredByString
     */
     void mapComplexs(GraphDatabaseService graphDb, String supportScoreString, String abundanceScoreString, String scoredByString, String mappedString){

        try (Transaction tx = graphDb.beginTx()) {

            ResourceIterator<Node> complexes = graphDb.findNodes(Label.label("Complex"));
            for (ResourceIterator<Node> it = complexes; it.hasNext(); ) {
                Node cplx = it.next();
                String mappedUIDs = "";
                HashSet<Node> components = new HashSet<>();
                HashSet<Node> nodes = recurseComponents(graphDb, cplx, components);
                String scoringUIDs = "";


                // for each node get score and avg it
                Double sScore = 0.0;
                Double eScore = 0.0;
                int size = nodes.size();

                // if one component is scored then score complex,
                // otherwise complex score should be left null
                boolean scoredComplex = false;
                for(Node component: nodes){
                    if(component.hasProperty(supportScoreString)){
                        scoredComplex = true;
                        Iterable<Relationship> compUIDRels = component.getRelationships(RelTypes.ID_BELONGS_TO);
                        for (Relationship compUIDRel: compUIDRels) {
                            Node uid = compUIDRel.getStartNode();
                            if (uid.hasProperty(mappedString)){
                                mappedUIDs = mappedUIDs.concat(uid.getProperty(PropertyType.UNIPROT_ID.toString()).toString()+ ", ");
                            }
                        }
                    }
                }

                if(scoredComplex){
                    for(Node component: nodes){
                        if (component.hasProperty(supportScoreString)){ // if its a scored prt
                            if(component.hasProperty(scoredByString)){ // multi UIDs
                                scoringUIDs = scoringUIDs + component.getProperty(scoredByString) + " , ";
                            }else{ // single UID
                                Iterable<Relationship> uidRELs = component.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                                if(getLength(uidRELs) > 1){
                                    for(Relationship uidRel: uidRELs){
                                        Node uidNode = uidRel.getStartNode();
                                        String uid = uidNode.getProperty(PropertyType.UNIPROT_ID.toString()).toString();

                                        scoringUIDs = scoringUIDs + uid + " , ";
                                    }
                                }
                            }
                            sScore +=  Double.parseDouble(component.getProperty(supportScoreString).toString());
                            eScore += Double.parseDouble(component.getProperty(abundanceScoreString).toString());
                        }
                    }
                    sScore = sScore/size; // TODO sscore is diluted by # of things in complex?
                    eScore = eScore/size; //TODO should escore be the avg of only things that have score?

                    cplx.setProperty(supportScoreString, (Math.round(sScore*1000d))/1000d);
                    cplx.setProperty(abundanceScoreString, (Math.round(eScore*1000d))/1000d);
                    cplx.setProperty(mappedString, mappedUIDs);
                }
            }

            tx.success();
        }

    }

    /**
     * Gets all base components of a complex (recurses complexes of complexes)
     * @param graphDb
     * @param complex
     * @param components
     * @return
     */
    private  HashSet<Node> recurseComponents(GraphDatabaseService graphDb, Node complex, HashSet<Node> components){
        // Components -> if the components are complexes or have members themselves then recurse
        if (complex.hasLabel(Label.label("Complex"))){
            //System.out.println("complex: " + entity + " components " + ((Complex) entity).getComponent());
            HashSet<Node> entityComponentSet = getComponent(complex);
            //components.addAll(entityComponentSet);
            for(Node physicalEntity: entityComponentSet){ // for each component of the current complex
                if(physicalEntity.hasLabel(Label.label("Complex"))){ // if it's component is a complex itself
                    recurseComponents(graphDb, physicalEntity,components); // recurse and get its components
                }else{
                    components.add(physicalEntity); // else add it's component to the list of components
                }
            }
        }
        return components;
    }

    /**
     * Gets all components of a complex (all components only one step away)
     * @param complex
     * @return
     */
    private HashSet<Node> getComponent( Node complex){
        HashSet<Node> components = new HashSet<>();
        Iterable<Relationship> relationships = complex.getRelationships(Direction.INCOMING, RelTypes.COMPONENT);
        for(Relationship rel: relationships){
            Node startNode = rel.getStartNode();
            components.add(startNode);
        }
        return components;
    }

    /**
     * adding weights to relationships
     */
    void addRelWeights(GraphDatabaseService graphDb, String supportScoreString, Double max, Double min){
        String tempSupportScoreString = supportScoreString.replaceAll("SUPPORT_SCORE_", "");
        String abundanceScoreString = "ABUNDANCE_SCORE_" + tempSupportScoreString;
        String weightSupportString = "WEIGHT_SUPPORT_" + tempSupportScoreString;
        String weightAbundanceString = "WEIGHT_ABUNDANCE_" + tempSupportScoreString;

        // for abundance weights take absolute value (will worrk for LFC and abundances)
        // for support weights, range is 0 - 1.5 and should just be what the support score is

        try (Transaction tx = graphDb.beginTx()) {

            // make abun.weight = the highest score in that experiment for all
            // make support weight = 1.5
            // because traversals minimize weight and we dont want these edges prioritized

            Double maxAbsScore = Math.max(Math.abs(min), Math.abs(max));
            ResourceIterable<Relationship> allRelationships = graphDb.getAllRelationships();
            for(Relationship rel: allRelationships){
                rel.setProperty(weightAbundanceString, (Math.round(maxAbsScore*1000d))/1000d);
                rel.setProperty(weightSupportString, 1.5);
            }


            // assign weights to incoming rlthnshps for pe nodes
            // once again weights cant be negative
            ResourceIterator<Node> pes = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            for (ResourceIterator<Node> it = pes; it.hasNext(); ) {
                Node pe = it.next();
                Iterable<Relationship> rltnshps = pe.getRelationships(Direction.INCOMING);
                for(Relationship rel: rltnshps){
                    if (pe.hasProperty(abundanceScoreString)){
                        double w = Double.parseDouble(pe.getProperty(abundanceScoreString).toString());
                        w = maxAbsScore - Math.abs(w);
                        rel.setProperty(weightAbundanceString, (Math.round(w*1000d))/1000d);

                        double ss = Double.parseDouble(pe.getProperty(supportScoreString).toString());
                        rel.setProperty(weightSupportString, (Math.round((1.5-ss)*1000d))/1000d);
                    }
                }
            }
            tx.success();
        }
    }

    /**
     * @throws IOException
     */
    public void printScoreDistributions() throws IOException{
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        //Create File to write outputs to
        FileWriter fstream  = new FileWriter(outputFile + "/scoreDistributions.tsv");
        BufferedWriter out = new BufferedWriter(fstream);
        out.write("UID\tDB_ID\tScore\teScore\tintegrated\tLocation\tPhosphorylated\tNumPhosns_pf\tPathway\n" );

        try (Transaction tx = graphDb.beginTx()) {
            ResourceIterator<Node> nodeSet = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));

            while(nodeSet.hasNext()){
                Node node = nodeSet.next();

                //get UID
                String UID = "";
                if(node.hasProperty(PropertyType.SCORED_BY.toString())){
                    UID = node.getProperty(PropertyType.SCORED_BY.toString()).toString();
                } else if (node.hasProperty(PropertyType.UNIPROT_ID.toString())){
                    UID = node.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                }else{
                    Iterable<Relationship> uidRels = node.getRelationships(RelTypes.ID_BELONGS_TO, Direction.INCOMING);
                    if(getLength(uidRels) == 1){
                        for(Relationship uidRel: uidRels){
                            UID = uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                        }
                    }else{
                        for(Relationship uidRel: uidRels){
                            UID = UID + uidRel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString() + ",";
                        }
                    }
                }

                //get Pathway
                String pathway = "-";
                for (Relationship pthRel : node.getRelationships(Direction.INCOMING, RelTypes.PATHWAY_COMPONENT)) {
                    pathway = pthRel.getStartNode().getProperty(PropertyType.DISPLAY_NAME.toString()).toString();
                }

                String SSCORE = "";
                String INTEGRATED = "";
                String LOCATION = "";
                String PHOSD = "";
                String NUM_PHOSNS = "";
                String  ESCORE = "";
                if(node.hasProperty(PropertyType.SUPPORT_SCORE.toString())){
                    SSCORE = node.getProperty(PropertyType.SUPPORT_SCORE.toString()).toString();
                    ESCORE = node.getProperty(PropertyType.ABUNDANCE_SCORE.toString()).toString();
                    if(node.hasProperty(PropertyType.INTEGRATED.toString())){
                        INTEGRATED = "1";
                    }else{
                        INTEGRATED = "0";
                    }
                    if(node.hasProperty(PropertyType.LOCATION.toString())){
                        LOCATION = node.getProperty(PropertyType.LOCATION.toString()).toString();
                    }else{
                        LOCATION = "NA";
                    }
                    Iterable<Relationship> phosRels = node.getRelationships(Direction.INCOMING, RelTypes.PHOSPHORYLATION);
                    if (getLength(phosRels) >0){
                        PHOSD = "1";
                    }else{
                        PHOSD = "0";
                    }
                    NUM_PHOSNS = getLength(phosRels).toString();
                    out.write(UID+ "\t" +
                            node.getId() + "\t" +
                            SSCORE +"\t" +
                            ESCORE +"\t" +
                            INTEGRATED +"\t" +
                            LOCATION +"\t" +
                            PHOSD + "\t" +
                            NUM_PHOSNS +"\t" +
                            pathway +"\n");
                }
            }
            tx.success();
        }
        out.close();
        graphDb.shutdown();
    }


    /**
     * @return
     */
    private Integer getLength(ResourceIterator<Object> thing) {
        Integer count = 0;
        for (ResourceIterator<Object> it = thing; it.hasNext(); ) {
            Object object = it.next();
            count ++;
        }
        return count;
    }

    /**
     * @param thing
     * @return
     */
    private Integer getLength(Iterable<Relationship> thing) {
        Integer count = 0;
        for (Relationship relationship : thing) {
            count++;
        }
        return count;
    }

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
