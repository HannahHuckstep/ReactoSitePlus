package org.wehi.hucksteph;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.biopax.paxtools.model.level2.protein;
import org.neo4j.graphalgo.GraphAlgoFactory;
import org.neo4j.graphalgo.PathFinder;
import org.neo4j.graphalgo.WeightedPath;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.graphdb.traversal.Evaluation;
import org.neo4j.graphdb.traversal.Evaluator;
import org.neo4j.graphdb.traversal.TraversalDescription;
import org.neo4j.graphdb.traversal.Traverser;
import org.neo4j.unsafe.impl.batchimport.input.InputException;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.Math.*;

/// THIS FILE NOT TO BE ADDED TO GIT
public class qPhosDatabase extends MeasuredDatabase{
    private String species;
    private final String UID_PATTERN = "[OPQ][0-9][A-Z0-9]{3}[0-9](\\-[0-9*]{1,2})?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(\\-[0-9*]{1,2})?";
    private final String humanUniProt = "https://www.uniprot.org/uniprot/?query=organism:9606&format=fasta&include=yes";
    private final String mouseUniProt = "https://www.uniprot.org/uniprot/?query=organism:10090&format=fasta&include=yes";

    public qPhosDatabase(File databaseDir, File outputFile) {
        super(databaseDir, outputFile);
    }

    public qPhosDatabase(File databaseDir) {
        super(databaseDir);
    }

    public void allDS() throws IOException {

        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        FileWriter fstream = new FileWriter(outputFile + "/TraversalReport_upstream.tsv");
        BufferedWriter out = new BufferedWriter(fstream);
        out.write("UniProtID\tUniProtNodeID\tProteoformName\tProteoformID\ttraversalDepth\n");

        FileWriter fstream1 = new FileWriter(outputFile + "/TraversalReport_downstream.tsv");
        BufferedWriter out1 = new BufferedWriter(fstream1);
        out1.write("UniProtID\tUniProtNodeID\tProteoformName\tProteoformID\ttraversalDepth\n");

        try (Transaction tx = graphDb.beginTx()) {
            ResourceIterator<Node> uids = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            ResourceIterator<Node> uids4len = graphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));

            TraversalDescription ds = graphDb.traversalDescription().breadthFirst()
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
                    .relationships(RelTypes.INPUT, Direction.OUTGOING)
                    .relationships(RelTypes.OUTPUT, Direction.OUTGOING)
                    .relationships(RelTypes.CONTROLS, Direction.OUTGOING)
                    .relationships(RelTypes.CATALYSIS, Direction.OUTGOING)
                    .relationships(RelTypes.COMPONENT, Direction.OUTGOING)
                    .relationships(RelationshipType.withName("ACTIVATION"), Direction.OUTGOING)
                    .relationships(RelationshipType.withName("INHIBITION"), Direction.OUTGOING);

            TraversalDescription us = graphDb.traversalDescription().breadthFirst()
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
                    .relationships(RelTypes.COMPONENT, Direction.INCOMING)
                    .relationships(RelationshipType.withName("ACTIVATION"), Direction.INCOMING)
                    .relationships(RelationshipType.withName("INHIBITION"), Direction.INCOMING);

            Integer count = 0;
            Integer i = getLength(uids4len);
            int onePercent = round(i/ 100);
            for (ResourceIterator<Node> it = uids; it.hasNext(); ) {
                Node uidNode = it.next();
                count ++;
//                if((count%onePercent) == 0){
//                    System.out.print("\rProgress: "+ (count/i)*100);
//                }

                Iterable<Relationship> relationships = uidNode.getRelationships(RelTypes.ID_BELONGS_TO);

                for (Relationship relationship : relationships) {
                    Node prot = relationship.getEndNode();

                    Traverser dsNodeTraverser = ds.traverse(prot);
                    Traverser usNodeTraverser = us.traverse(prot);

                    Integer dsDepth = 0;
                    for (Path nodePath : dsNodeTraverser) {

                        dsDepth++;
                    }

                    Integer usDepth = 0;
                    for (Path nodePath : usNodeTraverser) {
                        usDepth++;
                    }
                    //out.write("\tUniProtID\tUniProtNodeID\tProteoformName\tProeofromID\tDisplayName\ttraversalDepth\n");
                    out.write(uidNode.getProperty(PropertyType.DISPLAY_NAME.toString()) +"\t"+
                            uidNode.getId()  +"\t"+
                            prot.getProperty(PropertyType.DISPLAY_NAME.toString())   +"\t"+
                            prot.getId()   +"\t"+
                            (usDepth-1)  +"\n"
                    );

                    out1.write(uidNode.getProperty(PropertyType.DISPLAY_NAME.toString()) +"\t"+
                            uidNode.getId()  +"\t"+
                            prot.getProperty(PropertyType.DISPLAY_NAME.toString())   +"\t"+
                            prot.getId()   +"\t"+
                            (dsDepth-1)  +"\n"
                    );
                }
            }
            tx.success();
        }
        out.close();
        out1.close();
        graphDb.shutdown();

    }

    private Integer getLength(ResourceIterator<Node> thing) {
        Integer count = 0;
        while (thing.hasNext()){
            Node next = thing.next();
            count++;
        }
        return count;
    }

    // this one is to map all data points, outputs file to load inot proportion plot R script
    public void allMappings(File path2phosPeps) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        FileWriter fstream = new FileWriter(outputFile + "/qPhosAllMappings.tsv");
        BufferedWriter out = new BufferedWriter(fstream);
        out.write("Experiment\t" +
                "PepsMappedfromData\t" +
                "PeptidesNotMappedFromData\t" +
                "UIDsMappedFromData\t" +
                "UIDSNotMappedFromData\t" +
                "PhosnsMappedTo\t" +
                "PhosnsNotMappedTo\t" +
                "ProteoformsMappedTo\t" +
                "ProteoformsNotMappedTo\t" +
                "ComplexesMappedTo\t" +
                "ComplexesNotMappedTo\t" +
                "IntegratedMapped\t" +
                "IntegratedNotMapped\t" +
                "MaxAbundance\t" +
                "MinAbundance\n");


        String LRP_str = "UniProt_accession";
        String mod_pep = "Raw_peptide";
        String aVal = "Log2Ratio";
        String experiment = "Experiment";


        // read in human or mouse fasta file
        species = "human";
        HashMap<String, String> fastaMap = new HashMap<>();
        if (species.equalsIgnoreCase("human")) {
            // read in human fasta
            System.out.println("Reading in fasta");
            fastaMap = fasta2dict(humanUniProt);
        } else if (species.equalsIgnoreCase("mouse")) {
            // read in mouse fasta
            fastaMap = fasta2dict(mouseUniProt);
        }


        // Read in data into Dict
        BufferedReader BR = new BufferedReader(new FileReader(path2phosPeps));
        String line = "";
        int ID = 0;


        Integer leadingProtCol = 0;
        Integer modPepCol = 0;
        Integer aScoreCol = 0;
        Integer experimentCol = 0;

        Boolean lpBool = true;
        Boolean modPepBool = true;
        Boolean aSBool = true;
        Boolean expBool = true;

        // first get all experiemnt names by reading through the file
        HashSet<String> experiments = new HashSet<>();
        while ((line = BR.readLine()) != null) {
            ID++;
            String[] columns = line.split("\t");
            List<String> cols = new ArrayList<>();

            Integer colLen = columns.length;
            for (int i = 0; i < colLen; i++) {
                cols.add(columns[i]);
                if (columns[i].equals(experiment)) {
                    experimentCol = i;
                    expBool = false;
                } else if (columns[i].equalsIgnoreCase(LRP_str)) {
                    leadingProtCol = i;
                    lpBool = false;
                } else if (columns[i].equals(mod_pep)) {
                    modPepCol = i;
                    modPepBool = false;
                } else if (columns[i].equals(aVal)) {
                    aScoreCol = i;
                    aSBool = false;
                }
            }

            // Make sure the columns input exist
            if (lpBool) {
                throw new IOException("There is no " + LRP_str + " column in the input data file" +
                        "\nThese are the columns available \n" + cols);

            } else if (modPepBool) {
                throw new IOException("There is no " + mod_pep + " column in the input data file" +
                        "\nThese are the columns available \n" + cols);
            } else if (aSBool) {
                throw new IOException("There is no " + aVal + " column in the input data file" +
                        "\nThese are the columns available \n" + cols);
            } else if (expBool) {
                throw new IOException("There is no " + experiment + " column in the input data file" +
                        "\nThese are the columns available \n" + cols);
            }

            if (!columns[experimentCol].equals(experiment)) {
                experiments.add(columns[experimentCol]);
            }
        }
        System.out.println("Experiments to map: " + experiments.size());

        // next for each different experiment name, only gather the data in those experiments.
        int i = 0;
        for (String exprName : experiments) {
            i++;

            //Create File to write inputs/Outputs to
            System.out.println("\nExperiment #:" + i);
            System.out.println("Currently Mapping: " + exprName);
            String supportScoreString = "SUPPORT_SCORE_" + exprName;
            String abundanceScoreString = "ABUNDANCE_SCORE_" + exprName;
            String scoredByString = "SCORED_BY" + exprName;
            String mappedString = "MAPPED" + exprName;

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
            Integer pepNum = 0;


            // get columns of interest from data file
            while ((line = BR.readLine()) != null) {
                ID++;
                String[] columns = line.split("\t");

                if (columns[experimentCol].equals(exprName)) { // only per experiment
                    pepNum ++;
                    // extract LRP's
                    if (ID2LRP.containsKey(columns[leadingProtCol])) {
                        HashSet<String> temp = ID2LRP.get(columns[leadingProtCol]);
                        temp.add(String.valueOf(ID));
                        ID2LRP.put(columns[leadingProtCol], temp);
                    } else {
                        HashSet<String> temp = new HashSet<>();
                        temp.add(String.valueOf(ID));
                        ID2LRP.put(columns[leadingProtCol], temp);
                    }

                    //extract aScores:
                    if (columns[aScoreCol].isEmpty() | columns[aScoreCol].equalsIgnoreCase(aVal)) {
                        ID2eScore.put(String.valueOf(ID), null);
                    } else {
                        ID2eScore.put(String.valueOf(ID), Double.valueOf(columns[aScoreCol]));
                    }

                    if (!columns[modPepCol].equalsIgnoreCase(mod_pep)) {
                        // generate Start
                        String start = generateStartOfPep(fastaMap, columns[leadingProtCol], columns[modPepCol]);
                        if (start.isEmpty()) {
                            continue;
                        }

                        ID2Start.put(String.valueOf(ID), start);


                        // generate End
                        String end = generateEndOfPep(start, columns[modPepCol]);
                        ID2End.put(String.valueOf(ID), end);

                        // generate Top Mods
                        HashSet<String> mods = generateModList(ID, columns[modPepCol], ID2Start);
                        ID2TopMods.put(String.valueOf(ID), mods);
                    }
                }
            }

//            System.out.println("LRP: "+ ID2LRP);
//            System.out.println("ascore"+ ID2eScore);
//            System.out.println("id2start: " + ID2Start);
//            System.out.println("id2end: " + ID2End);
//            System.out.println("modlist: " + ID2TopMods);



            try (Transaction tx = graphDb.beginTx()) {

                Integer proteoformsScored = 0;
                double pepsSpanningModTemp = 0;
                double pepsSpanningMod = 0;
                Integer pepsMatched = 0;
                Integer UIDsMatched = 0;
                Integer UIDsNotMatched = 0;
                Integer PepsNotMatched = 0;

                // for each LRP in the data

                for (String LRP : ID2LRP.keySet()) {
                    HashSet<String> PeptideIDs = ID2LRP.get(LRP);

                    // get the UID from the DB
                    Node UID = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), LRP);


                    if (UID == null) {
                        UIDsNotMatched++;
                        PepsNotMatched += PeptideIDs.size();
                        continue;
                    } else {
                        UIDsMatched++;
                        pepsMatched += PeptideIDs.size();

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
                                if (!templocation.equals("unknown") & !templocation.equals("NA")) {
                                    Integer location = Integer.valueOf(phos.getProperty(PropertyType.LOCATION.toString()).toString());
                                    String type = phos.getProperty(PropertyType.TYPE.toString()).toString().replaceAll("p\\_", "");
                                    String rxmMod = type + "_" + location;
                                    UIDmodList.add(rxmMod);
                                    uidModList.put(location, type);
                                }
                            }
                        }


                        for (Relationship protRltnshp : protRelationships) {
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

                                if (!templocation.equals("unknown") & !templocation.equals("NA")) {
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
                                    protein.setProperty(supportScoreString, proteoformScore);
                                    protein.setProperty(scoredByString, LRP);
                                }
                            }else{ // no previous score - set it
                                protein.setProperty(supportScoreString, proteoformScore);
                            }

                            //System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()) +" pf sScore:"+ protein.getProperty(supportScoreString));

                            pepsSpanningModTemp = pepsSpanningModTemp/getLength(protRelationships);
                            pepsSpanningMod += pepsSpanningModTemp;
                            pepsSpanningModTemp =0;

                            //////////////////////////////////////////////////////////////////////////////////////////
                            // Setting abundance score:
                            // find the Peptide ID's with the biggest score

                            HashMap<String, Double> maxMin = new HashMap<>();
                            maxMin = abundanceScoreAvgHighestSupport(graphDb, pepID2sScore, ID2eScore, protein, abundanceScoreString, maxAbundance, minAbundance);

                            //System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()) +" pf aScore:"+ protein.getProperty(abundanceScoreString));

                            maxAbundance = maxMin.get("Max");
                            minAbundance = maxMin.get("Min");


                            //////////////////////////////////////////////////////////////////////////////////////////

                        }
                    }
                }

                out.write(exprName + "\t" + pepsMatched + "\t" + (pepNum - pepsMatched) + "\t" + UIDsMatched + "\t" + (UIDsNotMatched - 1) + "\t");

                tx.success();
            }
            mapComplexs(graphDb, supportScoreString, abundanceScoreString, scoredByString,mappedString);
            addRelWeights(graphDb, supportScoreString, maxAbundance, minAbundance);

            mappingReport(out, graphDb, supportScoreString, minAbundance, maxAbundance);

            resetScores(graphDb, exprName);
        }// expr for loop
        out.close();
        graphDb.shutdown();
    }

    //this one is for specifying an experiment to map -no output, used for MCN
    public void allMappings(File path2phosPeps, String expr) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);



        String LRP_str = "UniProt_accession";
        String mod_pep = "Raw_peptide";
        String aVal = "Log2Ratio";
        String experiment = "Experiment";


        // read in human or mouse fasta file
        species = "human";
        HashMap<String, String> fastaMap = new HashMap<>();
        if (species.equalsIgnoreCase("human")) {
            // read in human fasta
            System.out.println("Reading in fasta");
            fastaMap = fasta2dict(humanUniProt);
        } else if (species.equalsIgnoreCase("mouse")) {
            // read in mouse fasta
            fastaMap = fasta2dict(mouseUniProt);
        }


        // Read in data into Dict
        BufferedReader BR = new BufferedReader(new FileReader(path2phosPeps));
        String line = "";
        int ID = 0;


        Integer leadingProtCol = 0;
        Integer modPepCol = 0;
        Integer aScoreCol = 0;
        Integer experimentCol = 0;

        Boolean lpBool = true;
        Boolean modPepBool = true;
        Boolean aSBool = true;
        Boolean expBool = true;

        // first get all experiemnt names by reading through the file
        HashSet<String> experiments = new HashSet<>();
        while ((line = BR.readLine()) != null) {
            ID++;
            String[] columns = line.split("\t");
            List<String> cols = new ArrayList<>();

            Integer colLen = columns.length;
            for (int i = 0; i < colLen; i++) {
                cols.add(columns[i]);
                if (columns[i].equals(experiment)) {
                    experimentCol = i;
                    expBool = false;
                } else if (columns[i].equalsIgnoreCase(LRP_str)) {
                    leadingProtCol = i;
                    lpBool = false;
                } else if (columns[i].equals(mod_pep)) {
                    modPepCol = i;
                    modPepBool = false;
                } else if (columns[i].equals(aVal)) {
                    aScoreCol = i;
                    aSBool = false;
                }
            }

            // Make sure the columns input exist
            if (lpBool) {
                throw new IOException("There is no " + LRP_str + " column in the input data file" +
                        "\nThese are the columns available \n" + cols);

            } else if (modPepBool) {
                throw new IOException("There is no " + mod_pep + " column in the input data file" +
                        "\nThese are the columns available \n" + cols);
            } else if (aSBool) {
                throw new IOException("There is no " + aVal + " column in the input data file" +
                        "\nThese are the columns available \n" + cols);
            } else if (expBool) {
                throw new IOException("There is no " + experiment + " column in the input data file" +
                        "\nThese are the columns available \n" + cols);
            }

            if (!columns[experimentCol].equals(experiment)) { // if not first line
                if(columns[experimentCol].equals(expr)){    ///////////////////////////// ONLY MAP THIS EXPERIMENT
                    experiments.add(columns[experimentCol]);
                }
            }
        }
        System.out.println("Experiment Mapping  " + expr);

        // next for each different experiment name, only gather the data in those experiments.
        int i = 0;
        for (String exprName : experiments) {
            i++;

            //Create File to write inputs/Outputs to
            System.out.println("\nExperiment #:" + i);
            System.out.println("Currently Mapping: " + exprName);
            String supportScoreString = "SUPPORT_SCORE_" + exprName;
            String abundanceScoreString = "ABUNDANCE_SCORE_" + exprName;
            String scoredByString = "SCORED_BY" + exprName;
            String mappedString = "MAPPED" + exprName;

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
            Integer pepNum = 0;


            // get columns of interest from data file
            while ((line = BR.readLine()) != null) {
                ID++;
                String[] columns = line.split("\t");

                if (columns[experimentCol].equals(exprName)) { // only per experiment
                    pepNum ++;
                    // extract LRP's
                    if (ID2LRP.containsKey(columns[leadingProtCol])) {
                        HashSet<String> temp = ID2LRP.get(columns[leadingProtCol]);
                        temp.add(String.valueOf(ID));
                        ID2LRP.put(columns[leadingProtCol], temp);
                    } else {
                        HashSet<String> temp = new HashSet<>();
                        temp.add(String.valueOf(ID));
                        ID2LRP.put(columns[leadingProtCol], temp);
                    }

                    //extract aScores:
                    if (columns[aScoreCol].isEmpty() | columns[aScoreCol].equalsIgnoreCase(aVal)) {
                        ID2eScore.put(String.valueOf(ID), null);
                    } else {
                        ID2eScore.put(String.valueOf(ID), Double.valueOf(columns[aScoreCol]));
                    }

                    if (!columns[modPepCol].equalsIgnoreCase(mod_pep)) {
                        // generate Start
                        String start = generateStartOfPep(fastaMap, columns[leadingProtCol], columns[modPepCol]);
                        if (start.isEmpty()) {
                            continue;
                        }

                        ID2Start.put(String.valueOf(ID), start);


                        // generate End
                        String end = generateEndOfPep(start, columns[modPepCol]);
                        ID2End.put(String.valueOf(ID), end);

                        // generate Top Mods
                        HashSet<String> mods = generateModList(ID, columns[modPepCol], ID2Start);
                        ID2TopMods.put(String.valueOf(ID), mods);
                    }
                }
            }

//            System.out.println("LRP: "+ ID2LRP);
//            System.out.println("ascore"+ ID2eScore);
//            System.out.println("id2start: " + ID2Start);
//            System.out.println("id2end: " + ID2End);
//            System.out.println("modlist: " + ID2TopMods);



            try (Transaction tx = graphDb.beginTx()) {

                Integer proteoformsScored = 0;
                double pepsSpanningModTemp = 0;
                double pepsSpanningMod = 0;
                Integer pepsMatched = 0;
                Integer UIDsMatched = 0;
                Integer UIDsNotMatched = 0;
                Integer PepsNotMatched = 0;

                // for each LRP in the data

                for (String LRP : ID2LRP.keySet()) {
                    HashSet<String> PeptideIDs = ID2LRP.get(LRP);

                    // get the UID from the DB
                    Node UID = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), LRP);


                    if (UID == null) {
                        UIDsNotMatched++;
                        PepsNotMatched += PeptideIDs.size();
                        continue;
                    } else {
                        UIDsMatched++;
                        pepsMatched += PeptideIDs.size();

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
                                if (!templocation.equals("unknown") & !templocation.equals("NA")) {
                                    Integer location = Integer.valueOf(phos.getProperty(PropertyType.LOCATION.toString()).toString());
                                    String type = phos.getProperty(PropertyType.TYPE.toString()).toString().replaceAll("p\\_", "");
                                    String rxmMod = type + "_" + location;
                                    UIDmodList.add(rxmMod);
                                    uidModList.put(location, type);
                                }
                            }
                        }


                        for (Relationship protRltnshp : protRelationships) {
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

                                if (!templocation.equals("unknown") & !templocation.equals("NA")) {
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
                                    protein.setProperty(supportScoreString, proteoformScore);
                                    protein.setProperty(scoredByString, LRP);
                                }
                            }else{ // no previous score - set it
                                protein.setProperty(supportScoreString, proteoformScore);
                            }

                            //System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()) +" pf sScore:"+ protein.getProperty(supportScoreString));

                            pepsSpanningModTemp = pepsSpanningModTemp/getLength(protRelationships);
                            pepsSpanningMod += pepsSpanningModTemp;
                            pepsSpanningModTemp =0;

                            //////////////////////////////////////////////////////////////////////////////////////////
                            // Setting abundance score:
                            // find the Peptide ID's with the biggest score

                            HashMap<String, Double> maxMin = new HashMap<>();
                            maxMin = abundanceScoreAvgHighestSupport(graphDb, pepID2sScore, ID2eScore, protein, abundanceScoreString, maxAbundance, minAbundance);

                            //System.out.println(protein.getProperty(PropertyType.DISPLAY_NAME.toString()) +" pf aScore:"+ protein.getProperty(abundanceScoreString));

                            maxAbundance = maxMin.get("Max");
                            minAbundance = maxMin.get("Min");


                            //////////////////////////////////////////////////////////////////////////////////////////

                        }
                    }
                }



                tx.success();
            }
            mapComplexs(graphDb, supportScoreString, abundanceScoreString, scoredByString, mappedString);
            addRelWeights(graphDb, supportScoreString, maxAbundance, minAbundance);

        }// expr for loop
        System.out.println("Done Mapping" + expr);
        graphDb.shutdown();
    }

    private void mappingReport(BufferedWriter out2, GraphDatabaseService graphDb, String supportScoreString, double minAbundance, double maxAbundance) throws IOException {

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
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while (proteins.hasNext()){
                Node protein = proteins.next();
                if(protein.hasProperty(supportScoreString)){
                    pfsMappedTo++;
                }else{
                    pfsNotMappedTo++;
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

            out2.write( phosnsMappedTo+"\t" +
                    phosnsNotMappedTo+"\t" +
                    pfsMappedTo+"\t" +
                    pfsNotMappedTo+"\t" +
                    cplxsMappedTo+"\t" +
                    cplxsNotMappedTo+"\t" +
                    intsMappedTo+"\t" +
                    intsNotMappedTo+"\t" +
                    minAbundance+"\t" +
                    maxAbundance+"\n");
            tx.success();
        }
    }

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

        // peps in this format: KDsGLYLKELIEPVLtCFNDADsRL
        for(int i = 0; i <pep.length();i++ ){
            if(pep.charAt(i) == 's' | pep.charAt(i) == 't'| pep.charAt(i) == 'y'){
                char modType = pep.charAt((i)); // bc the char before is the phosphorylated one
                String x = ID2Start.get(String.valueOf(ID));
                int modLoc = (i-1) + Integer.valueOf(x);
                String modTypeStr = String.valueOf(modType).toUpperCase();
                String modLocStr = String.valueOf(modLoc);
                String mod = modTypeStr +"_"+ modLocStr;
                mods.add(mod);
            }
        }
        return mods;
    }

    private String generateStartOfPep(HashMap<String, String> fasta, String LRP, String modPep){
        String start = "";

        // get fasta seq
        String seq = "";
        if(fasta.containsKey(LRP)){
            seq = fasta.get(LRP);
        }

        // get mod pep
        modPep = modPep.toUpperCase();

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

    private String generateEndOfPep(String start, String modPep){
        String end = "";

        if(!start.isEmpty()){
            int strt = Integer.valueOf(start);
            int length = modPep.length();
            int nd = strt + length-1;
            end = String.valueOf(nd);
        }

        return end;
    }

    Integer getLength(Iterable<Relationship> thing) {
        Integer count = 0;
        for (Relationship relationship : thing) {
            count++;
        }
        return count;
    }

    public void allNbhds( Integer depth) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        // write new files
        FileWriter fstream = new FileWriter(outputFile + "/BinomialNeighbourhoods.tsv");
        BufferedWriter out = new BufferedWriter(fstream);
        out.write("Experiment\tnumSig\tNumBonfSig\n");



        try(Transaction tx = graphDb.beginTx()){
            //get all experiments in db
            ResourceIterable<String> allPropertyKeys = graphDb.getAllPropertyKeys();
            List<String> scores = new ArrayList<>();
            for (String property: allPropertyKeys) {
                if (property.startsWith("SUPPORT_SCORE_")){
                    scores.add(property);
                }
            }

            // get number of things measured in db (per experiment)
            HashMap<String, Integer> numMeasured = new HashMap<>();
            double numInDb = 0;
            for (String scoreString: scores) {
                String supportScoreStr = scoreString;
                numMeasured.put(scoreString,0);

                ResourceIterator<Node> allPEs = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
                while (allPEs.hasNext()) {
                    Node innerPE = allPEs.next();
                    if (innerPE.hasLabel(Label.label("Protein")) | innerPE.hasLabel(Label.label("Complex"))) {
                        if (innerPE.hasProperty(supportScoreStr)) {
                            Integer numMeasuredCount = numMeasured.get(supportScoreStr);
                            numMeasuredCount++;
                            numMeasured.put(supportScoreStr, numMeasuredCount);
                        }
                        numInDb++;
                    }
                }
            }

                // Find all neighbourhoods
            ResourceIterator<Node> PEs = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));

            //record the num sig in here
            HashMap<String, Integer> bonfSig = new HashMap<>();
            HashMap<String, Integer> pvalSig = new HashMap<>();

            Integer count = 0;
            for (ResourceIterator<Node> it = PEs; it.hasNext(); ) {
                count ++;

                Node pe = it.next();
                if(pe.hasLabel(Label.label("Protein")) | pe.hasLabel(Label.label("Complex"))) {

                    HashSet<Node> nbhd = neighbourTraversal(pe, depth, graphDb);
                    System.out.println("Neighbouthood Number:" + count);

                    // now get the # things measured and calc binom pval for that nbhd
                    Integer scoreCount = 0;
                    for (String scoreString: scores) {
                        scoreCount ++;
                        String experiment = scoreString.replaceAll("SUPPORT_SCORE_", "");
                        String supportScoreStr = scoreString;
                        bonfSig.put(experiment, 0);
                        pvalSig.put(experiment, 0);
                        Integer dbMeasured = numMeasured.get(supportScoreStr);


                        double nbhdSize = 0;
                        int nbhdMeasured = 0;
                        for(Node node:nbhd){
                            if (node.hasProperty(supportScoreStr)){
                                if(Double.parseDouble(node.getProperty(supportScoreStr).toString()) > 0){
                                    nbhdMeasured ++;
                                }
                                nbhdSize ++;
                            }
                        }

                        double probability = nbhdSize/numInDb;
                        Double pval = 0.0;
                        if(nbhdMeasured != 0){
                            BinomialTest bn = new BinomialTest();
                            pval = bn.binomialTest(dbMeasured, nbhdMeasured, probability, AlternativeHypothesis.GREATER_THAN);
                        }

                        Double bonfCorrected= pval * dbMeasured;

                        if(pval < 0.05){
                            Integer countSig = pvalSig.get(experiment);
                            countSig ++;
                            pvalSig.put(experiment, countSig);
                        }

                        if(bonfCorrected < 0.05 & nbhdSize > 0.0 ){
                            Integer countSig = pvalSig.get(experiment);
                            countSig ++;
                            bonfSig.put(experiment, countSig);
                        }
                    }
                }
            }
            Set<String> experiments = pvalSig.keySet();
            for(String experiment: experiments){
                Integer numPvalSig = pvalSig.get(experiment);
                Integer numBonfSig = bonfSig.get(experiment);
                out.write(experiment+"\t"+numPvalSig + "\t" + numBonfSig+ "\n");
            }

            tx.success();
        }
        out.close();
        graphDb.shutdown();
    }

    public void qPhosMinimalConnectionNetwork(String experiment) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputFile = getOutputFile();
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);

        try (Transaction tx = graphDb.beginTx()) {

            //for(String experiment: experimentStrs){
                System.out.println("\nMCN FOR "+ experiment);
                FileWriter fstream = new FileWriter(outputFile + "/MinimalConnectionNetworkReport_"+experiment+".tsv");
                BufferedWriter out = new BufferedWriter(fstream);
                out.write("pf\tmappedPf\tphos\tmappedPhos\tcplx\tmappedCplx\tkins\tmappedKins\tIntegrated\tmappedInt\tUIDs\tRXNs\n");

                String supportScoreStr = "SUPPORT_SCORE_" + experiment;
                String weightString = "WEIGHT_ABUNDANCE_" + experiment;

                HashSet<Node> mapped = new HashSet<>();
                ResourceIterator<Node> nodes = graphDb.findNodes(Label.label("Protein")); //TODO are we sure?
                while (nodes.hasNext()){
                    Node node = nodes.next();
                    if(node.hasProperty(supportScoreStr)){
                        mapped.add(node);
                    }
                }
                System.out.println("Mapped Proteins: "+ mapped.size());


                PathFinder<WeightedPath> dijkstraDS = GraphAlgoFactory.dijkstra(
                        PathExpanders.forTypesAndDirections(
                                RelTypes.INPUT, Direction.OUTGOING,
                                RelTypes.OUTPUT, Direction.OUTGOING,
                                RelTypes.CONTROLS, Direction.OUTGOING,
                                RelTypes.CATALYSIS, Direction.OUTGOING,
                                RelTypes.ID_BELONGS_TO, Direction.OUTGOING,
                                RelTypes.ID_BELONGS_TO, Direction.INCOMING,
                                RelTypes.COMPONENT, Direction.OUTGOING, //TODO ??
                                RelationshipType.withName("ACTIVATION"), Direction.OUTGOING,
                                RelationshipType.withName("INHIBITION"), Direction.OUTGOING),
                        weightString,
                        1
                );

                PathExpanderBuilder builder = PathExpanderBuilder.allTypes(Direction.OUTGOING)
                        .remove(RelationshipType.withName(PropertyType.SMALL_MOL_EDGE.toString()))
                        .remove(RelTypes.ID_BELONGS_TO);
                PathFinder<Path> finder = GraphAlgoFactory.shortestPath(builder.build(), 10);

                // make a hashset for all nodes and all relationships
                HashSet<Node> mcnNodes = new HashSet<>();
                HashSet<Relationship> mcnRels = new HashSet<>();


                // find paths between all
                for (Node start: mapped){

                    HashSet<Integer> uidsStart = new HashSet<Integer>();


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

                Integer numProteins = 0;
                Integer numMappedProteins = 0;
                Integer numPhosdProteins = 0;
                Integer numMappedPhosdProteins = 0;
                Integer numComplexes = 0;
                Integer numMappedComplexes = 0;
                Integer numRXNs = 0;
                Integer numKinases = 0;
                Integer numMappedKinases = 0;
                Integer numIntegrated = 0;
                Integer numMappedIntegrated = 0;
                HashSet<String> UIDs = new HashSet<>();
                for (Node node: mcnNodes){
                    if(node.hasLabel(Label.label("Protein"))){
                        numProteins++;
                        if(node.hasProperty(supportScoreStr)){
                            numMappedProteins++;
                        }
                        if (node.hasProperty(PropertyType.KINASE.toString())){
                            numKinases++;
                            if(node.hasProperty(supportScoreStr)){
                                numMappedKinases++;
                            }
                        }
                        if(node.hasProperty(PropertyType.INTEGRATED.toString())){
                            numIntegrated++;
                            if(node.hasProperty(supportScoreStr)){
                                numMappedIntegrated++;
                            }
                        }
                        if (node.hasRelationship(RelTypes.PHOSPHORYLATION)){
                            numPhosdProteins++;
                            if(node.hasProperty(supportScoreStr)){
                                numMappedPhosdProteins++;
                            }
                        }

                        Iterable<Relationship> uidRelationships = node.getRelationships(Direction.INCOMING, RelationshipType.withName(RelTypes.ID_BELONGS_TO.toString()));
                        for(Relationship rel: uidRelationships){
                            String s = rel.getStartNode().getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                            UIDs.add(s);
                        }


                    }else if(node.hasLabel(Label.label("Complex"))){
                        numComplexes++;
                        if(node.hasProperty(supportScoreStr)){
                            numMappedComplexes++;
                        }
                        if(node.hasProperty(PropertyType.INTEGRATED.toString())){
                            numIntegrated++;
                            if(node.hasProperty(supportScoreStr)){
                                numMappedIntegrated++;
                            }
                        }
                    }else if (node.hasLabel(Label.label("BiochemicalReaction"))){
                        numRXNs++;
                    }
                }

                out.write(numProteins +"\t"+ numMappedProteins +"\t"+
                        numPhosdProteins+"\t"+numMappedPhosdProteins+"\t"+
                        numComplexes+"\t"+numMappedComplexes+"\t"+
                        numKinases+"\t"+numMappedKinases+"\t"+
                        numIntegrated+"\t"+numMappedIntegrated+"\t"+
                        UIDs.size()+"\t"+numRXNs+"\n");

                out.close();
            //}
            tx.success();
        }
        System.out.println("MCN Done");
        graphDb.shutdown();
    }

}
