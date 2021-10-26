package org.wehi.hucksteph;

import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ReactomeDatabase extends EmbeddedNeo4jDatabase{


    public ReactomeDatabase( File databaseDir, File outputFile) {
        super( databaseDir, outputFile);
    }

    public ReactomeDatabase( File databaseDir) {
        super( databaseDir);
    }

    /**
     * Integrates PhosphositePlus with Reactome
     * @param PSPDatabaseDir the built PhosphositePlus K-S biopax database
     * @throws IOException
     */
    public void IntegratePSP(File PSPDatabaseDir, String species) throws IOException {
        File databaseDir = getDatabaseDir();
        File outputDir =  getOutputFile();

        if(!outputDir.exists()){
            outputDir.mkdir();
        }



        FileWriter fstream = new FileWriter(outputDir + "/IntegrationStats.txt");
        BufferedWriter out = new BufferedWriter(fstream);

        HashMap<String, String> uniProtNames = getUniProtNames(PSPDatabaseDir);

        // get phosphoSites and their Controllers from PSP
        // pspKS -> {targetUID:{targetPhosphosite:[controllingKinaseUIDs]}}
        HashMap<String, HashMap<String, HashSet<String>>> pspKS = getPspPhosAndControllers(PSPDatabaseDir, species);

        // first pass match phosphosites and
        GraphDatabaseService rxmGraphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);


        try (Transaction tx = rxmGraphDb.beginTx()) {

            ResourceIterator<Node> speciesNodes = rxmGraphDb.findNodes(Label.label("SPECIES"));

            String graphSpecies = "";
            while(speciesNodes.hasNext()){
                Node speciesNode = speciesNodes.next();
                String speciesStr = speciesNode.getProperty("Species").toString();
                if(speciesStr.equalsIgnoreCase("human")){
                    graphSpecies = "Human";
                }else{
                    graphSpecies = "Mouse";
                }
            }

            String specStr = "";
            if(species.equalsIgnoreCase("h") | species.equalsIgnoreCase("human") | species.equals("9606")){
                specStr = "Human";
            }else if(species.equalsIgnoreCase("m") | species.equalsIgnoreCase("mouse") | species.equals("10090")){
                specStr = "Mouse";
            }

            if (!graphSpecies.equals(specStr)){
                throw new Error("Species in Reactome database is not equal to the species specified");
            }



            // DONT INTEGRATE IF ALREADY INTEGRATED
            for (String allPropertyKey : rxmGraphDb.getAllPropertyKeys()) {
                if (allPropertyKey.equalsIgnoreCase(PropertyType.INTEGRATED.toString())){
                    System.out.println("Database has already been integrated");
                    System.exit(0);
                }
            }



            Integer UIDsnotinrxm = 0;
            Integer matchCount = 0;
            HashSet<String> uidsCanAddTo = new HashSet<>();
            Integer complexUIDswithPSPcontroller = 0;
            HashSet<Long> proteinUIDswithPSPcontroller = new HashSet<>();
            HashSet<Long> proteinUIDswithINCORRECTPSPcontroller = new HashSet<>();
            Integer PEUIDswithPSPcontroller = 0;
            Integer numBchmRxnsWithNoCTRL = 0;
            Integer numPhosdProtsNotAttavhedToARXN = 0;
            Integer numSinglyPhosdRXNsWithNoController = 0;
            Integer numPhosphositesAdded = 0;

            // getting original number of phsophosites
            ResourceIterator<Node> phosNum = rxmGraphDb.findNodes(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
            while(phosNum.hasNext()){
                Node x = phosNum.next();
                numPhosphositesAdded ++;
            }


            Set<String> pspUIDs = pspKS.keySet();
            // for each targetUID get the psite hashmap
            for (String pspUID : pspUIDs) {

                //pspControlledPhosns ->  {targetPhosphosite:[controllingKinaseUIDs]}
                HashMap<String, HashSet<String>> pspControlledPhosns = pspKS.get(pspUID);
                HashMap<String, HashSet<String>> completedPhosns = new HashMap<>();

                // get corresponding UID in rxm
                Node rxmUID = rxmGraphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), pspUID);

                if (rxmUID != null) {
                    uidsCanAddTo.add(pspUID);
                    // get all proteins
                    Iterable<Relationship> rxmUIDProteoforms = rxmUID.getRelationships();
                    for (Relationship rxmUIDProteoform : rxmUIDProteoforms) {
                        Node proteoform = rxmUIDProteoform.getEndNode();
                        Iterable<Relationship> rxmPhosns = proteoform.getRelationships(Direction.INCOMING, RelTypes.PHOSPHORYLATION);
                        // get all phosns for that proteoform
                        HashSet<String> rxmPhosnSet = new HashSet<>();
                        if (getLength(rxmPhosns) > 0) { //phos'd protein
                            for (Relationship phsonRL : rxmPhosns) {
                                String rxmtype = phsonRL.getStartNode().getProperty(PropertyType.TYPE.toString()).toString();
                                String rxmLoc = phsonRL.getStartNode().getProperty(PropertyType.LOCATION.toString()).toString();
                                rxmtype = rxmtype.substring(2, 3);
                                String rxmPhos = rxmtype + "-" + rxmLoc;
                                rxmPhosnSet.add(rxmPhos);
                            }

                            if(rxmPhosnSet.size() == 1){ // singly phosd
                                matchCount ++;

                                for(String currentPhos: rxmPhosnSet){
                                    if(pspControlledPhosns.containsKey(currentPhos)){ // if it's in psp set ( if we have a controller for it)

                                        // check if rxm prot is involved in a rxn
                                        Iterable<Relationship> rxnOutputs = proteoform.getRelationships(Direction.INCOMING, RelTypes.OUTPUT);
                                        if(getLength(rxnOutputs) >0){
                                            for(Relationship rxnOutput: rxnOutputs){
                                                Node rxn = rxnOutput.getStartNode();
                                                // controls
                                                /////////////// Phosphorylations that have a controller  /////////////////////
                                                Iterable<Relationship> rxnActivation = rxn.getRelationships(Direction.INCOMING,
                                                        RelationshipType.withName("ACTIVATION"),
                                                        RelationshipType.withName("INHIBITION"),
                                                        RelTypes.CATALYSIS);
                                                if(getLength(rxnActivation) >0 ){
                                                    for(Relationship act: rxnActivation){
                                                        Node catRXN = act.getStartNode();
                                                        Iterable<Relationship> ctrlRxns = catRXN.getRelationships(Direction.INCOMING, RelTypes.CONTROLS);
                                                        if(getLength(ctrlRxns) >0){
                                                            for(Relationship ctrlRxn: ctrlRxns){
                                                                Node controller = ctrlRxn.getStartNode();

                                                                // these are the kinases controlling phosns in psp
                                                                HashSet<String> pspControllerKinases = pspControlledPhosns.get(currentPhos);

                                                                // now you have w/e is controlling the phos rxn

                                                                String controllerType = controller.getProperty(PropertyType.TYPE.toString()).toString();
                                                                if (controllerType.equals("Complex")){ // has UID's as property
                                                                    String complexUIDs = controller.getProperty(PropertyType.UNIPROT_ID.toString()).toString();

                                                                    // check if complex has the uid from kinase controlling rxn in psp
                                                                    for (String pspControllerKinase: pspControllerKinases){
                                                                        if (complexUIDs.contains(pspControllerKinase)){
                                                                            complexUIDswithPSPcontroller ++;
                                                                            // DOESNT COUNT? CHECK INPUTS AND OUTPUTS & ADD NEW RXN AND CONTROL
                                                                        }
                                                                    }
                                                                }else if ( controllerType.equals("Protein")){ // UID attached node
                                                                    Iterable<Relationship> uidRLTNSHPs = controller.getRelationships(Direction.INCOMING, RelTypes.ID_BELONGS_TO);
                                                                    for(Relationship uidRLTNSHP : uidRLTNSHPs){
                                                                        Node UIDnode = uidRLTNSHP.getStartNode();
                                                                        String UIDstring = UIDnode.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                                                                        for (String pspControllerKinase: pspControllerKinases){
                                                                            if(pspControllerKinase.equals(UIDstring)){
                                                                                proteinUIDswithPSPcontroller.add(catRXN.getId());

                                                                                HashSet<String> temp;
                                                                                if(completedPhosns.containsKey(currentPhos)){
                                                                                    temp = completedPhosns.get(currentPhos);
                                                                                }else {
                                                                                    temp = new HashSet<>();
                                                                                }
                                                                                temp.add(pspControllerKinase);
                                                                                completedPhosns.put(currentPhos, temp);

                                                                            }else{ // protein is not controller expected
                                                                                proteinUIDswithINCORRECTPSPcontroller.add(catRXN.getId());

                                                                                //// THIS K-S PAIR ALREADY EXISTs
                                                                                Node beforeNode = null;
                                                                                Iterable<Relationship> rxnInputRels = rxn.getRelationships(Direction.INCOMING, RelTypes.INPUT);
                                                                                if(getLength(rxnInputRels) == 2){ // protein and ATP
                                                                                    for(Relationship rxnInputRel: rxnInputRels){ // get the protein
                                                                                        Node startNode = rxnInputRel.getStartNode();
                                                                                        if(startNode.getProperty(PropertyType.TYPE.toString()).equals("Protein")){
                                                                                            beforeNode = startNode;
                                                                                            makeBchmRxn(rxmGraphDb, beforeNode, proteoform, pspControllerKinase );
                                                                                        }
                                                                                    }
                                                                                }

                                                                                //// THIS K-S PAIR HAS BEEN ADDED
                                                                                HashSet<String> temp;
                                                                                if(completedPhosns.containsKey(currentPhos)){
                                                                                    temp = completedPhosns.get(currentPhos);
                                                                                }else {
                                                                                    temp = new HashSet<>();
                                                                                }
                                                                                temp.add(pspControllerKinase);
                                                                                completedPhosns.put(currentPhos, temp);
                                                                            }
                                                                        }
                                                                    }
                                                                }else if (controllerType.equals("PhysicalEntity")){ // MIGHT have UID's as property
                                                                    // DOESNT COUNT? CHECK INPUTS AND OUTPUTS & ADD NEW RXN AND CONTROL


                                                                    if(controller.hasProperty(PropertyType.SET.toString())){
                                                                        String peSet = controller.getProperty(PropertyType.SET.toString()).toString();
                                                                        String[] split = peSet.split(",");
                                                                        for (String member: split){
                                                                            // check if complex has the uid from kinase controlling rxn in psp
                                                                            for (String pspControllerKinase: pspControllerKinases){
                                                                                if (member.contains(pspControllerKinase)){
                                                                                    PEUIDswithPSPcontroller ++;
                                                                                    // DOESNT COUNT? CHECK INPUTS AND OUTPUTS & ADD NEW RXN AND CONTROL
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }// NO controller
                                                        else{
                                                            numSinglyPhosdRXNsWithNoController ++;
                                                        }
                                                    }
                                                }// NO CAT RXN
                                                else{
                                                    numBchmRxnsWithNoCTRL ++;
                                                }

                                            }
                                        }else{ // not in a reaction

                                            numPhosdProtsNotAttavhedToARXN ++;
                                            // check if proteoform with no mods (and no other uids) exists in the same CL
                                            HashSet<Integer> beforeNodeSet = new HashSet<>();
                                            Iterable<Relationship> rxmUIDProteoforms2 = rxmUID.getRelationships();
                                            for(Relationship rxmUIDProteoform2: rxmUIDProteoforms2) {
                                                Node proteoform2 = rxmUIDProteoform2.getEndNode();
                                                if(getLength(proteoform2.getRelationships(Direction.INCOMING, RelTypes.PHOSPHORYLATION)) == 0
                                                        & getLength(proteoform2.getRelationships(Direction.INCOMING, RelTypes.ID_BELONGS_TO)) == 1
                                                        & proteoform.getProperty(PropertyType.LOCATION.toString()).equals(proteoform2.getProperty(PropertyType.LOCATION.toString()))){
                                                    beforeNodeSet.add(( (int) proteoform2.getId()));
                                                }
                                            }

                                            // if only one, make it before node
                                            Node beforeProteinNode = null;
                                            if(beforeNodeSet.size() == 1){
                                                for(Integer beforeProteinID : beforeNodeSet) {
                                                    beforeProteinNode = rxmGraphDb.getNodeById(beforeProteinID);
                                                }
                                                // if multi before options, make a new one
                                            }else{
                                                beforeProteinNode = rxmGraphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
                                                beforeProteinNode.addLabel(Label.label("Protein"));
                                                beforeProteinNode.setProperty(PropertyType.LOCATION.toString(), proteoform.getProperty(PropertyType.LOCATION.toString()));
                                                beforeProteinNode.setProperty(PropertyType.DISPLAY_NAME.toString(), rxmUID.getProperty(PropertyType.UNIPROT_ID.toString()));
                                                beforeProteinNode.setProperty(PropertyType.TYPE.toString(), "Protein");
                                                beforeProteinNode.setProperty(PropertyType.INTEGRATED.toString(), "true");
                                                beforeProteinNode.setProperty(PropertyType.UNIPROT_ID.toString(), rxmUID.getProperty(PropertyType.UNIPROT_ID.toString()));
                                                if(proteoform.hasProperty(PropertyType.UNIPROT_NAME.toString())){
                                                    beforeProteinNode.setProperty(PropertyType.UNIPROT_NAME.toString(), proteoform.getProperty(PropertyType.UNIPROT_NAME.toString()));
                                                }

                                                // add before protein to UID
                                                Relationship relationshipTo = rxmUID.createRelationshipTo(beforeProteinNode, RelTypes.ID_BELONGS_TO);
                                                relationshipTo.setProperty(PropertyType.INTEGRATED.toString(), "true");
                                            }
                                            HashSet<String> pspControllerKinases = pspControlledPhosns.get(currentPhos);
                                            // create biochem rxn for each controller
                                            for (String pspController: pspControllerKinases){
                                                makeBchmRxn(rxmGraphDb, beforeProteinNode, proteoform, pspController);

                                                //// THIS K-S PAIR HAS BEEN ADDED
                                                HashSet<String> temp;
                                                if(completedPhosns.containsKey(currentPhos)){
                                                    temp = completedPhosns.get(currentPhos);
                                                }else {
                                                    temp = new HashSet<>();
                                                }
                                                temp.add(pspController);
                                                completedPhosns.put(currentPhos, temp);
                                            }
                                        }
                                    } // not in psp set (dont know controller)
                                }
                            } // multi phos'd
                        }
                    }
                    HashMap<String, HashSet<String>> toAdd = new HashMap<>();
                    for (String phosn: pspControlledPhosns.keySet()){
                        if(completedPhosns.containsKey(phosn)){
                            HashSet<String> completedKinases = completedPhosns.get(phosn);
                            HashSet<String> toAddKinases = pspControlledPhosns.get(phosn);
                            toAddKinases.removeAll(completedKinases);
                            if(!toAddKinases.isEmpty()){
                                toAdd.put(phosn, toAddKinases);
                            }
                        }else{
                            toAdd.put(phosn, pspControlledPhosns.get(phosn));
                        }
                    }
                    makeRxnNodes(rxmGraphDb, rxmUID, toAdd);

                }else{
                    UIDsnotinrxm ++;

                    // make before, after and rxn
                    Node newUIDnode = rxmGraphDb.createNode(Label.label(LabelTypes.UNIPROT_ID.toString()));
                    newUIDnode.setProperty(PropertyType.DISPLAY_NAME.toString(), pspUID);
                    newUIDnode.setProperty(PropertyType.STATUS.toString(), "Current");
                    newUIDnode.setProperty(PropertyType.TYPE.toString(), "ProteinID");
                    newUIDnode.setProperty(PropertyType.UNIPROT_ID.toString(), pspUID );
                    newUIDnode.setProperty(PropertyType.DB_CONNECTION.toString(), ("https://www.uniprot.org/uniprot/" + pspUID));
                    newUIDnode.setProperty(PropertyType.INTEGRATED.toString(), "true");

                    if(uniProtNames.containsKey(pspUID)){
                        String name = uniProtNames.get(pspUID);
                        newUIDnode.setProperty(PropertyType.UNIPROT_NAME.toString(), name);
                    }


                    makeRxnNodes(rxmGraphDb, newUIDnode, pspControlledPhosns);
                }
            }

            Integer newNumPhosphositesAdded = 0;
            ResourceIterator<Node> newPhosNum = rxmGraphDb.findNodes(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
            while(newPhosNum.hasNext()){
                Node x = newPhosNum.next();
                newNumPhosphositesAdded ++;
            }

            newNumPhosphositesAdded -= numPhosphositesAdded;

            System.out.println("The number of proteins(UIDs) that have phosphorylations in both RXM and PSP: "+uidsCanAddTo.size());
            System.out.println("The number of proteoforms with a single phosphosite in both RXM and PSP: "+ matchCount);
            System.out.println("Complex's with a component that is same controller as psp (doesnt matter, not the same rxn): " + complexUIDswithPSPcontroller);
            System.out.println("Proteins's that is same controller as psp: " + proteinUIDswithPSPcontroller.size());
            System.out.println("Proteins's that is DIFFERENT controller as psp: " + proteinUIDswithINCORRECTPSPcontroller.size());
            System.out.println("Physical Entity with a UID that is same controller as psp: " + PEUIDswithPSPcontroller);
            System.out.println("The number of Biochemical Reactions that have a phosphorylated protein as an output but no catalysis : " +  numSinglyPhosdRXNsWithNoController);
            System.out.println("The number of Biochemical Reactions that have a phosphorylated protein as an output but no controller : " + numBchmRxnsWithNoCTRL);
            System.out.println("The number of Phosphorylated proteoforms not involved in any biochemical Reactions: " + numPhosdProtsNotAttavhedToARXN);
            System.out.println("The number of UID's not in RXM (added): " + UIDsnotinrxm);
            System.out.println("The number of phosphosites not in RXM (added): " + newNumPhosphositesAdded);

            out.write("The number of proteins(UIDs) that have phosphorylations in both RXM and PSP: "+uidsCanAddTo.size() + "\n");
            out.write("The number of proteoforms with a single phosphosite in both RXM and PSP: "+ matchCount+ "\n");
            out.write("Complex's with a component that is same controller as psp (doesnt matter, not the same rxn): " + complexUIDswithPSPcontroller+ "\n");
            out.write("Proteins's that is same controller as psp: " + proteinUIDswithPSPcontroller.size()+ "\n");
            out.write("Proteins's that is DIFFERENT controller as psp: " + proteinUIDswithINCORRECTPSPcontroller.size()+ "\n");
            out.write("Physical Entity with a UID that is same controller as psp: " + PEUIDswithPSPcontroller+ "\n");
            out.write("The number of Biochemical Reactions that have a phosphorylated protein as an output but no catalysis : " +  numSinglyPhosdRXNsWithNoController+ "\n");
            out.write("The number of Biochemical Reactions that have a phosphorylated protein as an output but no controller : " + numBchmRxnsWithNoCTRL+ "\n");
            out.write("The number of Phosphorylated proteoforms not involved in any biochemical Reactions: " + numPhosdProtsNotAttavhedToARXN+ "\n");
            out.write("The number of UID's not in RXM (added): " + UIDsnotinrxm+ "\n");
            out.write("The number of phosphosites not in RXM (added): " + newNumPhosphositesAdded+ "\n");


            tx.success();
            out.close();
        }
        rxmGraphDb.shutdown();

    }

    /**
     * Makes all nodes needed to form a reaction, input nodes, output nodes, phosphorylation nodes,
     * calls a function to make biochemical reaction nodes, and attaches nodes to UIDs.
     * @param rxmGraphDb
     * @param UIDnode
     * @param pspControlledPhosns
     */
    private void makeRxnNodes(GraphDatabaseService rxmGraphDb, Node UIDnode, HashMap<String, HashSet<String>> pspControlledPhosns){
        String pspUID = UIDnode.getProperty(PropertyType.DISPLAY_NAME.toString()).toString();

        Node beforeProteinNode = rxmGraphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
        beforeProteinNode.addLabel(Label.label("Protein"));
        beforeProteinNode.setProperty(PropertyType.LOCATION.toString(), "unknown");
        beforeProteinNode.setProperty(PropertyType.DISPLAY_NAME.toString(), pspUID);
        beforeProteinNode.setProperty(PropertyType.UNIPROT_ID.toString(), pspUID);
        beforeProteinNode.setProperty(PropertyType.TYPE.toString(), "Protein");
        beforeProteinNode.setProperty(PropertyType.INTEGRATED.toString(), "true");
        if (UIDnode.hasProperty(PropertyType.UNIPROT_NAME.toString())){
            beforeProteinNode.setProperty(PropertyType.UNIPROT_NAME.toString(), UIDnode.getProperty(PropertyType.UNIPROT_NAME.toString()));
        }


        Relationship relationshipTo = UIDnode.createRelationshipTo(beforeProteinNode, RelTypes.ID_BELONGS_TO);
        relationshipTo.setProperty(PropertyType.INTEGRATED.toString(), "true");


        for(String pspPhosn: pspControlledPhosns.keySet()){

            Node afterNode = rxmGraphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            afterNode.addLabel(Label.label("Protein"));
            afterNode.setProperty(PropertyType.LOCATION.toString(), "unknown");
            afterNode.setProperty(PropertyType.DISPLAY_NAME.toString(), pspUID +"_p-"+pspPhosn);
            afterNode.setProperty(PropertyType.UNIPROT_ID.toString(), pspUID);
            afterNode.setProperty(PropertyType.TYPE.toString(), "Protein");
            afterNode.setProperty(PropertyType.INTEGRATED.toString(), "true");
            if (UIDnode.hasProperty(PropertyType.UNIPROT_NAME.toString())){
                beforeProteinNode.setProperty(PropertyType.UNIPROT_NAME.toString(), UIDnode.getProperty(PropertyType.UNIPROT_NAME.toString()));
            }

            Relationship relationshipTo1 = UIDnode.createRelationshipTo(afterNode, RelTypes.ID_BELONGS_TO);
            relationshipTo1.setProperty(PropertyType.INTEGRATED.toString(), "true");


            Pattern locationPattern = Pattern.compile("-(.*)");
            Matcher mLoc = locationPattern.matcher(pspPhosn);
            mLoc.find();
            String location = mLoc.group(1);



            Pattern typePattern = Pattern.compile("S|T|Y");
            Matcher mType = typePattern.matcher(pspPhosn);
            mType.find();
            String type = mType.group(0);


            Node rxmPhosnOnAfterNode = rxmGraphDb.createNode(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
            rxmPhosnOnAfterNode.setProperty(PropertyType.DISPLAY_NAME.toString(), pspUID +"p_"+pspPhosn);
            rxmPhosnOnAfterNode.setProperty(PropertyType.LOCATION.toString(),location);
            rxmPhosnOnAfterNode.setProperty(PropertyType.TYPE.toString(), "p_"+type);
            rxmPhosnOnAfterNode.setProperty(PropertyType.DB_ID.toString(), pspUID +"_"+pspPhosn );
            rxmPhosnOnAfterNode.setProperty(PropertyType.INTEGRATED.toString(), "true");

            Relationship relationshipTo2 = rxmPhosnOnAfterNode.createRelationshipTo(afterNode, RelTypes.PHOSPHORYLATION);
            relationshipTo2.setProperty(PropertyType.INTEGRATED.toString(), "true");


            HashSet<String> pspControllerKinases = pspControlledPhosns.get(pspPhosn);

            for (String controller: pspControllerKinases){
                makeBchmRxn(rxmGraphDb, beforeProteinNode, afterNode, controller);
            }

        }
    }

    /**
     * Makes the biochemical reactions nodes and controllers when integrating PSP with Reactome
     * @param rxmGraphDb
     * @param beforeProtein
     * @param afterProtein
     * @param controller
     */
    private void makeBchmRxn( GraphDatabaseService rxmGraphDb ,Node beforeProtein, Node afterProtein, String controller){

        // make RXN,
        Node rxnNode = rxmGraphDb.createNode(Label.label(LabelTypes.INTERACTION.toString()));
        rxnNode.addLabel(Label.label("BiochemicalReaction"));
        rxnNode.setProperty(PropertyType.TYPE.toString(), "BiochemicalReaction");
        rxnNode.setProperty( PropertyType.DISPLAY_NAME.toString(), controller + " can Phosphorylate " +
                beforeProtein.getProperty(PropertyType.DISPLAY_NAME.toString()));
        rxnNode.setProperty(PropertyType.INTEGRATED.toString(), "true");


        Relationship relationshipTo2 = beforeProtein.createRelationshipTo(rxnNode, RelTypes.INPUT);
        relationshipTo2.setProperty(PropertyType.INTEGRATED.toString(), "true");
        Relationship relationshipTo3 = rxnNode.createRelationshipTo(afterProtein, RelTypes.OUTPUT);
        relationshipTo3.setProperty(PropertyType.INTEGRATED.toString(), "true");



        // make controller rxn Node
        Node catNode = rxmGraphDb.createNode(Label.label(LabelTypes.INTERACTION.toString()));
        catNode.addLabel(Label.label("Catalysis"));
        catNode.setProperty(PropertyType.DISPLAY_NAME.toString(), "Catalysis");
        catNode.setProperty(PropertyType.TYPE.toString(), "Catalysis");
        catNode.setProperty(PropertyType.INTEGRATED.toString(), "true");


        Relationship relationshipTo4 = catNode.createRelationshipTo(rxnNode, RelTypes.CATALYSIS);
        relationshipTo4.setProperty(PropertyType.INTEGRATED.toString(), "true");


        String location = afterProtein.getProperty(PropertyType.LOCATION.toString()).toString();

        // add ATP & ADP
        HashMap< String, Long> atpIDnLoc = new HashMap<>();
        HashMap< String, Long> adpIDnLoc = new HashMap<>();
        ResourceIterator<Node> atp = rxmGraphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()), PropertyType.DISPLAY_NAME.toString(), "ATP");
        ResourceIterator<Node> adp = rxmGraphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()), PropertyType.DISPLAY_NAME.toString(), "ADP");

        while(atp.hasNext()) {
            Node ATPnext = atp.next();
            for (Map.Entry<String, Object> stringObjectEntry : ATPnext.getProperties(PropertyType.LOCATION.toString()).entrySet()) {
                atpIDnLoc.put( stringObjectEntry.getValue().toString(), ATPnext.getId());
            }
        }
        while(adp.hasNext()) {
            Node ADPnext = adp.next();
            for (Map.Entry<String, Object> stringObjectEntry : ADPnext.getProperties(PropertyType.LOCATION.toString()).entrySet()) {
                adpIDnLoc.put( stringObjectEntry.getValue().toString(), ADPnext.getId());
            }
        }


        if((atpIDnLoc.containsKey(location) & adpIDnLoc.containsKey(location))  ){
            Node atpNode = rxmGraphDb.getNodeById(atpIDnLoc.get(location));
            Node adpNode = rxmGraphDb.getNodeById(adpIDnLoc.get(location));

            Relationship relationshipTo = atpNode.createRelationshipTo(rxnNode, RelTypes.INPUT);
            relationshipTo.setProperty(PropertyType.SMALL_MOL_EDGE.toString(), "Input");
            relationshipTo.setProperty(PropertyType.INTEGRATED.toString(), "true");

            Relationship relationshipTo1 = rxnNode.createRelationshipTo(adpNode, RelTypes.OUTPUT);
            relationshipTo1.setProperty(PropertyType.SMALL_MOL_EDGE.toString(), "Output");
            relationshipTo1.setProperty(PropertyType.INTEGRATED.toString(), "true");

        }else{ // no atp/adp in that location


            Node atpNode = rxmGraphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            atpNode.addLabel(Label.label("SmallMolecule"));
            atpNode.setProperty(PropertyType.DISPLAY_NAME.toString(), "ATP");
            atpNode.setProperty(PropertyType.LOCATION.toString(), "unknown");
            atpNode.setProperty(PropertyType.TYPE.toString(), "SmallMolecule");
            atpNode.setProperty(PropertyType.INTEGRATED.toString(), "true");


            Node adpNode = rxmGraphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            adpNode.addLabel(Label.label("SmallMolecule"));
            adpNode.setProperty(PropertyType.DISPLAY_NAME.toString(), "ADP");
            adpNode.setProperty(PropertyType.LOCATION.toString(), "unknown");
            adpNode.setProperty(PropertyType.TYPE.toString(), "SmallMolecule");
            adpNode.setProperty(PropertyType.INTEGRATED.toString(), "true");


            Relationship relationshipTo = atpNode.createRelationshipTo(rxnNode, RelTypes.INPUT);
            relationshipTo.setProperty(PropertyType.SMALL_MOL_EDGE.toString(), "Input");
            relationshipTo.setProperty(PropertyType.INTEGRATED.toString(), "true");

            Relationship relationshipTo1 = rxnNode.createRelationshipTo(adpNode, RelTypes.OUTPUT);
            relationshipTo1.setProperty(PropertyType.SMALL_MOL_EDGE.toString(), "Output");
            relationshipTo1.setProperty(PropertyType.INTEGRATED.toString(), "true");
        }


        // find or create a controller
        Node controllerNode = findOrCreateController( rxmGraphDb, controller, location);

        Relationship relationshipTo = controllerNode.createRelationshipTo(catNode, RelTypes.CONTROLS);
        relationshipTo.setProperty(PropertyType.INTEGRATED.toString(), "true");


    }

    /**
     * When reaction has a controller, either find an appropriate existing controller or create a new one
     * @param rxmGraphDb
     * @param UID
     * @param location
     * @return
     */
    private Node findOrCreateController( GraphDatabaseService rxmGraphDb, String UID, String location){

        Node controllerNode = null;

        Node UIDnode = rxmGraphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), UID);

        if(UIDnode != null){

            HashSet<Integer> potentialControllers = new HashSet<>();
            Iterable<Relationship> UIDrelationships = UIDnode.getRelationships();
            for(Relationship uidRelationship: UIDrelationships){
                Node proteoform = uidRelationship.getEndNode();

                // if controller is unmodified and in the same location as outgoing protein
                // add to list of potential controllers
                if(getLength(proteoform.getRelationships(Direction.INCOMING, RelTypes.PHOSPHORYLATION))== 0 &
                        proteoform.getProperty(PropertyType.LOCATION.toString()).equals(location)){
                    potentialControllers.add((int) proteoform.getId());
                }
            }

            if (potentialControllers.size() == 1){
                //  the perfect controller, no mods in right CL
                for(Integer perfectController: potentialControllers){
                    controllerNode = rxmGraphDb.getNodeById(perfectController);
                }
            }else{
                // no controllers with no mods in right CL
                // multi potential controllers with no mods in right CL
                // make new one
                // problem where we might make a new ctrl node for each rxn if we have more than one to start with
                controllerNode = rxmGraphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
                controllerNode.addLabel(Label.label("Protein"));
                controllerNode.setProperty(PropertyType.LOCATION.toString(), location);
                controllerNode.setProperty(PropertyType.TYPE.toString(), "Protein");
                controllerNode.setProperty(PropertyType.DISPLAY_NAME.toString(), UID);
                controllerNode.setProperty(PropertyType.UNIPROT_ID.toString(), UID);
                controllerNode.setProperty(PropertyType.INTEGRATED.toString(), "true");

                Relationship relationshipTo = UIDnode.createRelationshipTo(controllerNode, RelTypes.ID_BELONGS_TO);
                relationshipTo.setProperty(PropertyType.INTEGRATED.toString(), "true");


            }

            // UID for psp Controller doesnt exist
        }else{
            // make a new UID and protein for it in right CL
            Node newUIDnode = rxmGraphDb.createNode(Label.label(LabelTypes.UNIPROT_ID.toString()));
            newUIDnode.setProperty(PropertyType.DISPLAY_NAME.toString(), UID);
            newUIDnode.setProperty(PropertyType.STATUS.toString(), "Current");
            newUIDnode.setProperty(PropertyType.TYPE.toString(), "ProteinID");
            newUIDnode.setProperty(PropertyType.UNIPROT_ID.toString(), UID );
            newUIDnode.setProperty(PropertyType.DB_CONNECTION.toString(), ("https://www.uniprot.org/uniprot/" + UID));
            newUIDnode.setProperty(PropertyType.INTEGRATED.toString(), "true");


            controllerNode = rxmGraphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            controllerNode.addLabel(Label.label("Protein"));
            controllerNode.setProperty(PropertyType.LOCATION.toString(), location);
            controllerNode.setProperty(PropertyType.TYPE.toString(), "Protein");
            controllerNode.setProperty(PropertyType.DISPLAY_NAME.toString(), UID);
            controllerNode.setProperty(PropertyType.UNIPROT_ID.toString(), UID);
            controllerNode.setProperty(PropertyType.INTEGRATED.toString(), "true");


            Relationship relationshipTo = newUIDnode.createRelationshipTo(controllerNode, RelTypes.ID_BELONGS_TO);
            relationshipTo.setProperty(PropertyType.INTEGRATED.toString(), "true");

        }

        return controllerNode;
    }

    /**
     * creates a nested hashmap containing all of PSP's phosphorylations and their respective controllers
     * @param PSPDatabaseDir
     * @return
     */
    private HashMap<String,HashMap<String, HashSet<String>>> getPspPhosAndControllers(File PSPDatabaseDir, String species){
        GraphDatabaseService pspGraphDb = new GraphDatabaseFactory().newEmbeddedDatabase(PSPDatabaseDir);

        // pspKS -> {targetUID:{targetPhosphosite:[controllingKinaseUIDs]}}
        HashMap<String,HashMap<String, HashSet<String>>> pspKS = new HashMap<>();
        try(Transaction tx = pspGraphDb.beginTx()){


            ResourceIterator<Node> speciesNodes = pspGraphDb.findNodes(Label.label("SPECIES"));

            String graphSpecies = "";
            while(speciesNodes.hasNext()){
                Node speciesNode = speciesNodes.next();
                String speciesStr = speciesNode.getProperty("Species").toString();
                if(speciesStr.equalsIgnoreCase("human")){
                    graphSpecies = "Human";
                }else{
                    graphSpecies = "Mouse";
                }
            }

            String specStr = "";
            if(species.equalsIgnoreCase("h") | species.equalsIgnoreCase("human") | species.equals("9606")){
                specStr = "Human";
            }else if(species.equalsIgnoreCase("m") | species.equalsIgnoreCase("mouse") | species.equals("10090")){
                specStr = "Mouse";
            }

            if (!graphSpecies.equals(specStr)){
                throw new Error("Species in PhosphoSitePlus database is not equal to the species specified");
            }


            // get all UIDs
            ResourceIterator<Node> pspUIDs = pspGraphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            // for each UID find all psites
            while (pspUIDs.hasNext()){
                Node pspUID = pspUIDs.next();
                if(pspUID.getProperty(PropertyType.STATUS.toString()).toString().equalsIgnoreCase("current")) {

                    String currentUID = pspUID.getProperty(PropertyType.UNIPROT_ID.toString()).toString();

                    // get all proteins for that uid
                    for (Relationship uidRLTNSHP : pspUID.getRelationships()) {
                        Node target = uidRLTNSHP.getEndNode();

                        // if protein is phosphorylated
                        Iterable<Relationship> targetPhosRLTNSHPs = target.getRelationships(Direction.INCOMING, RelTypes.PHOSPHORYLATION);
                        Integer length = getLength(targetPhosRLTNSHPs);

                        if (length > 0) {
                            String currentPhos = "";
                            for (Relationship targetPhosRLTNSHP : targetPhosRLTNSHPs) {
                                String type = targetPhosRLTNSHP.getStartNode().getProperty(PropertyType.TYPE.toString()).toString();
                                type = type.substring(2, 3);
                                String loc = targetPhosRLTNSHP.getStartNode().getProperty(PropertyType.LOCATION.toString()).toString();
                                currentPhos = type + "-" + loc;
                            }
                            // have to traverse all the way up to kinase UID
                            Iterable<Relationship> rxnRLTNSHPs = target.getRelationships(Direction.INCOMING, RelTypes.OUTPUT);
                            for (Relationship rxnRLTNSHP : rxnRLTNSHPs) {
                                Iterable<Relationship> catRLTNSHPs = rxnRLTNSHP.getStartNode().getRelationships(Direction.INCOMING, RelTypes.CATALYSIS);
                                for (Relationship catRLTNSHP : catRLTNSHPs) {
                                    Iterable<Relationship> ctrlRLTNSHPs = catRLTNSHP.getStartNode().getRelationships(Direction.INCOMING, RelTypes.CONTROLS);
                                    for (Relationship ctrlRLTNSHP : ctrlRLTNSHPs) {
                                        Iterable<Relationship> kinaseUIDRLTNSHPs = ctrlRLTNSHP.getStartNode().getRelationships(Direction.INCOMING, RelTypes.ID_BELONGS_TO);
                                        for (Relationship kinaseUIDRLTNSHP : kinaseUIDRLTNSHPs) {
                                            String controllingKinaseUID = kinaseUIDRLTNSHP.getStartNode().getProperty(PropertyType.DISPLAY_NAME.toString()).toString();

                                            //System.out.println(currentUID + " " + currentPhos + " controlled by" + controllingKinaseUID);
                                            HashMap<String, HashSet<String>> temp = new HashMap<>();// phosns: [contrilling kinases]
                                            if (pspKS.containsKey(currentUID)){
                                                temp = pspKS.get(currentUID);
                                                HashSet<String> temp2 = new HashSet<>(); // [controlling kinases]
                                                if(temp.containsKey(currentPhos)){
                                                    temp2 = temp.get(currentPhos);
                                                }
                                                temp2.add(controllingKinaseUID);
                                                temp.put(currentPhos,temp2);
                                            }else{
                                                HashSet<String> temp2 = new HashSet<>(); // [controlling kinases]
                                                temp2.add(controllingKinaseUID);
                                                temp.put(currentPhos,temp2);
                                            }
                                            pspKS.put(currentUID, temp);
                                        }
                                    }

                                }
                            }
                        } //else this uid belongs to a non phos'd protein (master or kinase/phosphatase)
                    }
                } //else UID NOT correct species

            }

            tx.success();
        }
        pspGraphDb.shutdown();
        return pspKS;
    }

    /**
     * Gets the UniProt names from Reactome, or PSP and adds them to integrated nodes
     * @param PSPDatabaseDir
     * @return
     */
    private HashMap<String, String> getUniProtNames(File PSPDatabaseDir){
        GraphDatabaseService pspGraphDb = new GraphDatabaseFactory().newEmbeddedDatabase(PSPDatabaseDir);
        HashMap<String, String> IDs2Names = new HashMap<>();

        try(Transaction tx = pspGraphDb.beginTx()){

            ResourceIterator<Node> pspsUIDs = pspGraphDb.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            while(pspsUIDs.hasNext()){
                Node pspUID = pspsUIDs.next();
                if(pspUID.hasProperty(PropertyType.UNIPROT_NAME.toString())){
                    String ID = pspUID.getProperty(PropertyType.UNIPROT_ID.toString()).toString();
                    String name = pspUID.getProperty(PropertyType.UNIPROT_NAME.toString()).toString();
                    IDs2Names.put(ID, name);
                }
            }
            tx.success();
        }
        pspGraphDb.shutdown();
        return  IDs2Names;
    }

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

