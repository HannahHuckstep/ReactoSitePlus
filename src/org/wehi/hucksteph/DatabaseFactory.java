package org.wehi.hucksteph;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.time.StopWatch;
import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level2.conversion;
import org.biopax.paxtools.model.level2.entity;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Entity;
import org.biopax.paxtools.model.level3.Process;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.unsafe.impl.batchimport.input.InputException;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DatabaseFactory {
    private File owlFile;
    private File databaseDir;
    private boolean update;
    private String species;

    private  Node createdNode;
    private Node current;
    private GraphDatabaseService graphDb;
    private String xmlBase;

    private final String secondaryAccessionFile = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/sec_ac.txt";
    private final String humanUniProtNames = "https://www.uniprot.org/uniprot/?query=*&format=tab&columns=id,genes(PREFERRED)&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes";
    private final String humanUIDs = "https://www.uniprot.org/uniprot/?query=organism:9606&format=tab&columns=id";
    private final String mouseUniProtNames = "https://www.uniprot.org/uniprot/?query=*&format=tab&columns=id,genes(PREFERRED)&fil=organism:%22Mus%20musculus%20(Mouse)%20[10090]%22%20AND%20reviewed:yes";
    private final String mouseUIDs = "https://www.uniprot.org/uniprot/?query=organism:10090&format=tab&columns=id";
    private final String UID_PATTERN_NO_ISO = "[OPQ][0-9][A-Z0-9]{3}[0-9](\\-[0-9*]{1,2})?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9])";
    private final String UID_PATTERN = "[OPQ][0-9][A-Z0-9]{3}[0-9](\\-[0-9*]{1,2})?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(\\-[0-9*]{1,2})?";
    private final String KINASES = "https://www.uniprot.org/uniprot/?query=keyword:%22Kinase%20[KW-0418]%22&format=list&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22";
    private final String TRANSCRIPTION_FACTORS = "https://www.uniprot.org/uniprot/?query=goa:(%22DNA-binding%20transcription%20factor%20activity%20[3700]%22)%20(reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22)&format=list";
    private final String CELL_SURFACE_RECEPTORS = "https://www.uniprot.org/uniprot/?query=goa:(%22cell%20surface%20receptor%20signaling%20pathway%20involved%20in%20cell-cell%20signaling%20[1905114]%22)%20(reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22)&format=list";
    private final String KINASES_MOUSE = "https://www.uniprot.org/uniprot/?query=keyword:%22Kinase%20[KW-0418]%22&format=list&fil=organism:%22Mus%20musculus%20(Mouse)%20[10090]%22";
    private final String TRANSCRIPTION_FACTORS_MOUSE = "https://www.uniprot.org/uniprot/?query=goa:(%22DNA-binding%20transcription%20factor%20activity%20[3700]%22)%20(reviewed:yes%20organism:%22Mus%20musculus%20(Mouse)%20[10090]%22)&format=list";
    private final String CELL_SURFACE_RECEPTORS_MOUSE = "https://www.uniprot.org/uniprot/?query=goa:(%22cell%20surface%20receptor%20signaling%20pathway%20involved%20in%20cell-cell%20signaling%20[1905114]%22)%20(reviewed:yes%20organism:%22Mus%20musculus%20(Mouse)%20[10090]%22)&format=list";




    public DatabaseFactory(File owlFile, File databaseDir, boolean update, String species) {
        if (!FilenameUtils.getExtension(owlFile.toString()).equalsIgnoreCase("owl")) {
            throw new IllegalArgumentException("input file not an owl file");
        } else {
            this.owlFile = owlFile;
            this.databaseDir = databaseDir;
            this.update = update;
            this.species = species;

            if (databaseDir.exists()) {
                Scanner scanner = new Scanner(System.in);

                System.out.println("*** WARNING ***");
                System.out.println("The folder "+ databaseDir+" already exists. \n"+
                        "To ensure a new database is made, this program will delete this folder: " + databaseDir +
                        "\nand generate a new one with the same name.\n" +
                        "Would you like to proceed?(y/n)");

                // delete existing database
                if (scanner.next().equalsIgnoreCase("y")) {
                    try {
                        FileUtils.deleteDirectory(databaseDir);
                    } catch (IOException e) {
                        e.printStackTrace();
                        System.exit(1);
                    }
                    // create new folder and database
                    databaseDir.mkdir();
                    graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
                    registerShutdownHook(graphDb);

                } else{
                    System.out.println("Database creation aborted");
                    System.exit(1);
                }

            } else {
                // folder did not already exist so make new folder and database
                databaseDir.mkdir();
                graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDir);
                registerShutdownHook(graphDb);


            }
        }
    }

    public DatabaseFactory(File owlFile, File databaseDir) {
        this(owlFile, databaseDir,true,"h");
    }

    public DatabaseFactory(String speciesStr) {
        this.species = speciesStr;
    }

    public DatabaseFactory() {
    }

    /**
     *  createDBfromOWL creates a neo4j embedded graph database
     *  updates uniprot id's automatically using the latest version of uniprot if updates parameter set to true
     */
    public void createDBfromOWL(){

        databaseDir.mkdir();

        // Create Model to get all conversions
        FileInputStream fin = null;
        try {
            fin = new FileInputStream(owlFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(1);
        }
        BioPAXIOHandler handler = new SimpleIOHandler();
        Model model = handler.convertFromOWL(fin);
        xmlBase = model.getXmlBase();

        if(update){
            model = updateBiopaxModel(model);
            xmlBase = model.getXmlBase();
        }

        Set<Conversion> conversionSet = model.getObjects(Conversion.class);
        Set<TemplateReaction> templateReactionSet = model.getObjects(TemplateReaction.class);
        Set<TemplateReactionRegulation> templateReactionRegulationSet = model.getObjects(TemplateReactionRegulation.class);
        Set<Protein> proteinSet = model.getObjects(Protein.class);
        Set<Pathway> pathwaySet = model.getObjects(Pathway.class);

        try (Transaction tx = graphDb.beginTx()) {
            // set species node
            if(species.equalsIgnoreCase("h") | species.equalsIgnoreCase("human") | species.equals("9606")){
                Node species = graphDb.createNode(Label.label("SPECIES"));
                species.setProperty("Species", "Human");

            }else if(species.equalsIgnoreCase("m") | species.equalsIgnoreCase("mouse") | species.equals("10090")){
                Node species = graphDb.createNode(Label.label("SPECIES"));
                species.setProperty("Species", "Mouse");
            }
            else{
                throw new InputException("Species must be either Human (\"h\", \"human\", \"9606\"), or Mouse (\"m\", \"mouse\", \"10090\")");
            }
            tx.success();
        }


        try (Transaction tx = graphDb.beginTx()) {
            for (Conversion conversion : conversionSet) {
                // current reaction we are looking at
                Node rxn = createNode(conversion);

                //get all input, outputs and controls
                Set<PhysicalEntity> LeftSet = conversion.getLeft();
                Set<PhysicalEntity> RightSet = conversion.getRight();
                Set<Control> controllerSet = conversion.getControlledOf();

                // create all nodes if they don't already exist
                // create separate edge type for small molecules for traversals
                for (PhysicalEntity left : LeftSet) {
                    Node input = createNode(left);
                    if(input.getProperty(PropertyType.TYPE.toString()).toString().equals("SmallMolecule")){
                        //input.createRelationshipTo(rxn, RelTypes.SMALL_MOL_EDGE);
                        Relationship relationship = input.createRelationshipTo(rxn, RelTypes.INPUT);
                        relationship.setProperty(PropertyType.SMALL_MOL_EDGE.toString(), "Input");
                    }else{
                        input.createRelationshipTo(rxn, RelTypes.INPUT);
                    }
                }
                for (PhysicalEntity right : RightSet) {
                    Node output = createNode(right);
                    if(output.getProperty(PropertyType.TYPE.toString()).toString().equals("SmallMolecule")){
                        //rxn.createRelationshipTo(output, RelTypes.SMALL_MOL_EDGE);
                        Relationship relationship = rxn.createRelationshipTo(output, RelTypes.OUTPUT);
                        relationship.setProperty(PropertyType.SMALL_MOL_EDGE.toString(), "Output");
                    }else {
                        rxn.createRelationshipTo(output, RelTypes.OUTPUT);
                    }
                }
                for (Control controlRXN : controllerSet) {
                    Node control = createNode(controlRXN);
                    if (controlRXN.getControlType() == null){
                        control.createRelationshipTo(rxn, RelTypes.CATALYSIS);
                    }else {
                        control.createRelationshipTo(rxn, RelationshipType.withName(controlRXN.getControlType().toString()));
                    }
                    Set<Controller> controllersOfControl = controlRXN.getController();
                    for (Controller controllers : controllersOfControl) {
                        Node controller = createNode(controllers);
                        controller.createRelationshipTo(control, RelTypes.CONTROLS);

                    }
                }
            }
            tx.success();
        }


        try (Transaction tx = graphDb.beginTx()) {
            // add all proteins not in a reaction
            for(Protein protein: proteinSet) {
                Node protNode = graphDb.findNode(Label.label("Protein"), PropertyType.DB_ID.toString(), " ");
                if(protNode == null){
                    createNode(protein);
                }
            }
            // add template reactions and their controllers

            //controller -- ctrls --> templateReactionRegulation -- ACTIVATES --> template Reaction -- OUTPUT --> product
            for(TemplateReaction tr: templateReactionSet){
                Node trNode = graphDb.createNode(Label.label(LabelTypes.INTERACTION.toString()));
                trNode.addLabel(Label.label(classToString(tr.getClass().toString())));
                trNode.setProperty(PropertyType.DB_ID.toString(), entityToString(tr));
                if (tr.getName().isEmpty()) {
                    trNode.setProperty(PropertyType.DISPLAY_NAME.toString(), "NA");
                } else {
                    trNode.setProperty(PropertyType.DISPLAY_NAME.toString(), tr.getName().toString());
                }
                trNode.setProperty(PropertyType.TYPE.toString(), classToString(tr.getClass().toString()));
                trNode.setProperty(PropertyType.DB_CONNECTION.toString(), databaseLink(tr));

                Set<PhysicalEntity> products = tr.getProduct();
                for (PhysicalEntity product: products) {
                    Node productNode = graphDb.findNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()), PropertyType.DB_ID.toString(), entityToString(product));
                    trNode.createRelationshipTo(productNode, RelTypes.OUTPUT);
                }

            }

            for(TemplateReactionRegulation trr: templateReactionRegulationSet){
                Node trrNode = graphDb.createNode(Label.label(LabelTypes.INTERACTION.toString()));
                trrNode.addLabel(Label.label(classToString(trr.getClass().toString())));
                trrNode.setProperty(PropertyType.DB_ID.toString(), entityToString(trr));
                if (trr.getName().isEmpty()) {
                    trrNode.setProperty(PropertyType.DISPLAY_NAME.toString(), "NA");
                } else {
                    trrNode.setProperty(PropertyType.DISPLAY_NAME.toString(), trr.getName().toString());
                }
                trrNode.setProperty(PropertyType.TYPE.toString(), classToString(trr.getClass().toString()));
                trrNode.setProperty(PropertyType.DB_CONNECTION.toString(), databaseLink(trr));

                Set<Controller> controllerSet = trr.getController();
                for (Controller controller: controllerSet) {
                    Node ctrlNode = createNode(controller);
                    ctrlNode.createRelationshipTo(trrNode, RelTypes.CONTROLS);
                }
                Set<Process> controlledSet = trr.getControlled();
                ControlType controlType = trr.getControlType();
                for (Process rxn: controlledSet) {
                    Node rxnNode = graphDb.findNode(Label.label(LabelTypes.INTERACTION.toString()), PropertyType.DB_ID.toString(), entityToString(rxn));
                    trrNode.createRelationshipTo(rxnNode, RelationshipType.withName(controlType.toString()));
                }
            }
            tx.success();
        }


        try (Transaction tx = graphDb.beginTx()) {

            // make a node for each pathway
            for (Pathway pathway : pathwaySet) {
                Set<Xref> xrefs = pathway.getXref();
                String id= "";
                    for(Xref xref: xrefs){
                        if(xref instanceof Xref){
                            String db = xref.getDb();
                            if(db != null){
                                if (db.equals("Reactome")){
                                    id = (xref).getId();
                                }
                            }
                        }
                    }

                Node pth = graphDb.createNode(Label.label(LabelTypes.PATHWAY.toString()));
                pth.setProperty(PropertyType.DISPLAY_NAME.toString(), pathway.getDisplayName());
                pth.setProperty(PropertyType.DB_ID.toString(), id);
                pth.setProperty(PropertyType.DB_CONNECTION.toString(), databaseLink(pathway));
                pth.setProperty(PropertyType.TYPE.toString(), "Pathway");
            }

            // then go through and add attachments and components
            // adding connecting inputs, outputs and controllers to it
            for (Pathway pathway : pathwaySet) {
                // get path unique DB_ID
                Set<Xref> xrefs = pathway.getXref();
                String id= "";
                for(Xref xref: xrefs){
                    if(xref instanceof Xref){
                        String db = xref.getDb();
                        if(db != null){
                            if (db.equals("Reactome")){
                                id = (xref).getId();
                            }
                        }
                    }
                }

                // find the node
                Node pth = graphDb.findNode(Label.label(LabelTypes.PATHWAY.toString()), PropertyType.DB_ID.toString(), id);
                Set<Process> pathwayComponent = pathway.getPathwayComponent();

                // get interactions that are a part of each pathway
                for (Process process : pathwayComponent) {
                    if(process instanceof Interaction){
                        Node intraction = graphDb.findNode(Label.label(LabelTypes.INTERACTION.toString()),
                                PropertyType.DB_ID.toString(), entityToString(process));

                        // if the interaction node exists
                        if(!(intraction == null)){
                            if(!relationshipBetween(pth, intraction)){
                                pth.createRelationshipTo(intraction, RelTypes.PATHWAY_COMPONENT);
                            }

                            Iterable<Relationship> relationships = intraction.getRelationships();

                            // for each reaction get its inputs, outputs and controllers and connect them to the path node
                            for(Relationship relationship: relationships){
                                if(relationship.isType(RelTypes.INPUT)){
                                    Node inStartNode = relationship.getStartNode();
                                    if(!relationshipBetween(pth, inStartNode)){
                                        pth.createRelationshipTo(inStartNode, RelTypes.PATHWAY_COMPONENT);
                                    }
                                }else if(relationship.isType(RelTypes.OUTPUT)){
                                    Node outEndNode = relationship.getEndNode();
                                    if(!relationshipBetween(pth, outEndNode)){
                                        pth.createRelationshipTo(outEndNode, RelTypes.PATHWAY_COMPONENT);
                                    }
                                }else if(relationship.isType(RelationshipType.withName("ACTIVATION")) ||
                                        relationship.isType(RelationshipType.withName("INHIBITION"))){
                                    Node ctrlStartNode = relationship.getStartNode();
                                    Iterable<Relationship> ctrlStartNodeRelationships = ctrlStartNode.getRelationships(Direction.INCOMING);
                                    for(Relationship relationship1: ctrlStartNodeRelationships){
                                        if(relationship1.isType(RelTypes.CONTROLS)){
                                            Node controller = relationship1.getStartNode();
                                            if(!relationshipBetween(pth,controller)){
                                                pth.createRelationshipTo(controller, RelTypes.PATHWAY_COMPONENT);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if(process instanceof Pathway){

                        Set<Xref> xrefs2 = process.getXref();
                        String id2= "";
                        for(Xref xref: xrefs2){
                            if(xref instanceof Xref){
                                String db = xref.getDb();
                                if(db != null){
                                    if (db.equals("Reactome")){
                                        id2 = (xref).getId();
                                    }
                                }
                            }
                        }
                        Node componentPth = graphDb.findNode(Label.label(LabelTypes.PATHWAY.toString()), PropertyType.DB_ID.toString(), id2);
                        pth.createRelationshipTo(componentPth, RelTypes.SUB_PATHWAY);
                    }
                }
            }
            tx.success();
        }


        // adds physical entity components not directly involved in interactions
        try (Transaction tx = graphDb.beginTx()) {
            addComponents( model);
            if(species.equalsIgnoreCase("Human") || species.equalsIgnoreCase("h")){
                addLabelsToDB("KINASE", KINASES);
                addLabelsToDB("TRANSCRIPTION_FACTOR", TRANSCRIPTION_FACTORS);
                addLabelsToDB("CELL_SURFACE_RECEPTOR", CELL_SURFACE_RECEPTORS);
            }else if(species.equalsIgnoreCase("Mouse") | species.equalsIgnoreCase("m")){
                addLabelsToDB("KINASE", KINASES_MOUSE);
                addLabelsToDB("TRANSCRIPTION_FACTOR", TRANSCRIPTION_FACTORS_MOUSE);
                addLabelsToDB("CELL_SURFACE_RECEPTOR", CELL_SURFACE_RECEPTORS_MOUSE);
            }
            addUniProtGeneNames();
            tx.success();
        }

        graphDb.shutdown();

    }

    /**
     * Checks to see what kind of object entity is and if it already exists as a node in the database.
     * If it doesnt exist, the appropriate function to create that node.
     * @param entity
     * is the biopax entity object to be made into a node
     */
    private  Node createNode(Entity entity) {

        // if its a physical entity
        if (entity instanceof PhysicalEntity) {
            // if it's a protein make the node and a node for its uniprot ID & its phosphorylations
            if (entity instanceof Protein) {
                if (graphDb.findNode(Label.label(classToString(classToString(entity.getClass().toString()))), PropertyType.DB_ID.toString(), entityToString(entity)) == null) {
                    createdNode = createProtein(entity);
                } else {
                    createdNode = graphDb.findNode(Label.label(classToString(entity.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(entity));
                }
            } else {
                if (graphDb.findNode(Label.label(classToString(entity.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(entity)) == null) {
                    createdNode = createPENode(entity);
                } else {
                    createdNode = graphDb.findNode(Label.label(classToString(entity.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(entity));
                }
            }
        }

        // If it's an interaction
        else if (entity instanceof Interaction) {

            if (graphDb.findNode(Label.label(classToString(entity.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(entity) ) == null) {
                current = graphDb.createNode(Label.label(classToString(entity.getClass().toString())));
                current.addLabel(Label.label(LabelTypes.INTERACTION.toString()));
                current.setProperty(PropertyType.DB_ID.toString(), entityToString(entity));
                if (entity instanceof Catalysis) {
                    current.setProperty(PropertyType.DISPLAY_NAME.toString(), "Catalysis");
                } else if (entity.getName().isEmpty()) {
                    current.setProperty(PropertyType.DISPLAY_NAME.toString(), "NA");
                } else {
                    current.setProperty(PropertyType.DISPLAY_NAME.toString(), entity.getName().toString());
                }
                String typeString = classToString(entity.getClass().toString());
                current.setProperty(PropertyType.TYPE.toString(), typeString);

                String dbconnecton1 = databaseLink(entity);
                current.setProperty(PropertyType.DB_CONNECTION.toString(), dbconnecton1);

                createdNode = current;
            } else {
                createdNode = graphDb.findNode(Label.label(classToString(entity.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(entity));
            }

        }

        return createdNode;

    }

    /**
     * Takes in a database and a model (from CreateDB) and adds all (not interacting) components of a complex to the database
     * Called at the end of createDB (b/c createDB only adds objects involved in an interaction to the database)
     * @param model
     * the current biopax file
     */
    private  void addComponents( Model model) {

        try (Transaction tx = graphDb.beginTx()) {

            // get all complexes in graphDb
            ResourceIterator<Node> complexes = graphDb.findNodes(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()), PropertyType.TYPE.toString(), "Complex");

            // traverse them
            while (complexes.hasNext()) {

                Node complexDB = complexes.next();
                //System.out.println("Current Complex: "+ complexDB.getProperty(PropertyType.DB_ID.toString()));

                // get the corresponding complex in Reactome via the unique Reactome ID
                String dbid = model.getXmlBase() + complexDB.getProperty(PropertyType.DB_ID.toString()).toString();
                BioPAXElement complexRXMmodel = model.getByID(dbid);
                Complex complexRXM = (Complex) complexRXMmodel;

                // get all of the Complex's components
                Set<PhysicalEntity> componentSet = complexRXM.getComponent();
                for(PhysicalEntity physicalEntity: componentSet){
                    // if that component doesnt exist
                    if (graphDb.findNode(Label.label(classToString(physicalEntity.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(physicalEntity)) == null) {
                        if(physicalEntity instanceof Protein){
                            Node protComp = createProtein(physicalEntity);
                            //System.out.println("Component: " + protComp.getProperty(PropertyType.DISPLAY_NAME.toString()));
                            protComp.createRelationshipTo(complexDB, RelTypes.COMPONENT);

                        }else if(physicalEntity instanceof Complex){
                            Node complexComp = createPENode(physicalEntity);
                            //System.out.println("Component being made: " + complexComp.getProperty(PropertyType.DB_ID.toString()));
                            complexComp.createRelationshipTo(complexDB, RelTypes.COMPONENT);
                            if(!((Complex) physicalEntity).getComponent().isEmpty()){
                                recurseRXMComponents((Complex) physicalEntity);
                            }
                        } else{
                            Node peComp = createPENode(physicalEntity);
                            //System.out.println("Component: " + peComp.getProperty(PropertyType.DISPLAY_NAME.toString()));
                            peComp.createRelationshipTo(complexDB, RelTypes.COMPONENT);
                        }

                    }else{
                        Node component = graphDb.findNode(Label.label(classToString(physicalEntity.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(physicalEntity));
                        //System.out.println("Component exists: " + component.getProperty(PropertyType.DB_ID.toString()));
                        component.createRelationshipTo(complexDB, RelTypes.COMPONENT);
                    }
                }
                Set<PhysicalEntity> memberPhysicalEntitySet = complexRXM.getMemberPhysicalEntity();
                for(PhysicalEntity physicalEntity: memberPhysicalEntitySet){
                    if (graphDb.findNode(Label.label(classToString(physicalEntity.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(physicalEntity)) == null) {
                        if(physicalEntity instanceof Protein){
                            Node protComp = createProtein(physicalEntity);
                            //System.out.println("MemPhysEnt: " + protComp.getProperty(PropertyType.DISPLAY_NAME.toString()));
                            protComp.createRelationshipTo(complexDB, RelTypes.COMPONENT);

                        }else {
                            Node peComp = createPENode(physicalEntity);
                            //System.out.println("MemPhysEnt: " + peComp.getProperty(PropertyType.DISPLAY_NAME.toString()));
                            peComp.createRelationshipTo(complexDB, RelTypes.COMPONENT);

                        }
                    }else{
                        Node component = graphDb.findNode(Label.label(classToString(physicalEntity.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(physicalEntity));
                        //System.out.println("MemPhysEnt: " + component.getProperty(PropertyType.DISPLAY_NAME.toString()));
                        component.createRelationshipTo(complexDB, RelTypes.COMPONENT);
                    }
                }
            }

            //TODO get all PE's and traverse and add members
            tx.success();
        }
    }

    /**
     * if an entity is a complex, recurse its components and add it to the database
     * @param complex
     * a biopax complex entity object
     */
    private  void recurseRXMComponents(Complex complex){

        // get the current complex node
        Node currentComplex = graphDb.findNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()), PropertyType.DB_ID.toString(), entityToString(complex));
        //System.out.println("\tCurrent Recurse Complex: " + currentComplex.getProperty(PropertyType.DB_ID.toString()));
        // get the list of components
        Set<PhysicalEntity> componentList = complex.getComponent();

        // for each one
        for (PhysicalEntity component: componentList) {
            // if the node doesn't exist make it
            if (graphDb.findNode(Label.label(classToString(component.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(component)) == null) {

                if (component instanceof Protein) {
                    Node protComp = createProtein(component);
                    //System.out.println("\tRecuse making Component: " + protComp.getProperty(PropertyType.DB_ID.toString()));
                    protComp.createRelationshipTo(currentComplex, RelTypes.COMPONENT);

                    // if the component is a complex itself recurse down the list of sub-components and attach them all to their super-complex
                } else if (component instanceof Complex ) {
                    Node complexComp = createPENode(component);
                    //System.out.println("\tRecurse making Component: " + complexComp.getProperty(PropertyType.DB_ID.toString()));
                    complexComp.createRelationshipTo(currentComplex, RelTypes.COMPONENT);
                    if (!((Complex) component).getComponent().isEmpty()){
                        recurseRXMComponents((Complex) component);
                    }
                } else {
                    Node peComp = createPENode(component);
                    //System.out.println("\tRecurse Component: " + peComp.getProperty(PropertyType.DB_ID.toString()));
                    peComp.createRelationshipTo(currentComplex, RelTypes.COMPONENT);

                }
                // the node does already exist
            }else{
                Node inDBNode = graphDb.findNode(Label.label(classToString(component.getClass().toString())), PropertyType.DB_ID.toString(), entityToString(component));
                //System.out.println("\tRecurse existing Component: " + inDBNode.getProperty(PropertyType.DB_ID.toString()));
                inDBNode.createRelationshipTo(currentComplex, RelTypes.COMPONENT);
            }
        }
    }

    /**
     * Creates a protein node in the database
     * @param entity
     * the model entity to be created
     * @return Node - the node to add to the database
     */
    private Node createProtein( Entity entity){

        current = graphDb.createNode(Label.label(classToString(entity.getClass().toString())));

        current.addLabel(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
        // unique ID
        current.setProperty(PropertyType.DB_ID.toString(), entityToString(entity));
        if (entity.getDisplayName() == null) {
            current.setProperty(PropertyType.DISPLAY_NAME.toString(), "NA");
        } else {
            current.setProperty(PropertyType.DISPLAY_NAME.toString(), entity.getDisplayName());
        }

        String typeString = classToString(entity.getClass().toString());
        current.setProperty(PropertyType.TYPE.toString(), typeString);

        if (((PhysicalEntity) entity).getCellularLocation() == null) {
            current.setProperty(PropertyType.LOCATION.toString(), "NA");
        } else {
            current.setProperty(PropertyType.LOCATION.toString(), entityCellularLocation(entity));
        }

        String dbconnecton1 = databaseLink(entity);
        current.setProperty(PropertyType.DB_CONNECTION.toString(), dbconnecton1);


        // if it is a protein and not a PE
        if (((Protein) entity).getMemberPhysicalEntity().isEmpty()) {
            // get UID(s)
            Set<String> uids = getUniprotID((Protein) entity);
            // if a there are UIDs
            creatUIDNode(uids);


            // get all modifications for a protein
            // get mods will check for mods on it's entity ref and add them to this prot
            ArrayList<String> mods = getMod((Protein) entity);
            if (!mods.isEmpty()) {
                // iterate through them
                for (String mod : mods) {
                    Pattern pattern = Pattern.compile("(?<=\\_).*$");
                    Matcher matcher = pattern.matcher(mod);
                    matcher.find();
                    // Location
                    if (mod.contains("p-S")) {
                        Node node = graphDb.createNode(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
                        node.setProperty(PropertyType.DB_ID.toString(), entityToString(entity) + "_" + mod);
                        node.setProperty(PropertyType.LOCATION.toString(), matcher.group(0));
                        node.setProperty(PropertyType.TYPE.toString(), "p_S");
                        node.setProperty(PropertyType.DISPLAY_NAME.toString(),  "p-S_"+ matcher.group(0));
                        node.createRelationshipTo(current, RelTypes.PHOSPHORYLATION);
                    } else if (mod.contains("p-T")) {
                        Node node = graphDb.createNode(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
                        node.setProperty(PropertyType.DB_ID.toString(), entityToString(entity) + "_" + mod);
                        node.setProperty(PropertyType.LOCATION.toString(), matcher.group(0));
                        node.setProperty(PropertyType.TYPE.toString(), "p_T");
                        node.setProperty(PropertyType.DISPLAY_NAME.toString(),  "p-T_"+ matcher.group(0));
                        node.createRelationshipTo(current, RelTypes.PHOSPHORYLATION);
                    } else if (mod.contains("p-Y")) {
                        Node node = graphDb.createNode(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
                        node.setProperty(PropertyType.DB_ID.toString(), entityToString(entity) + "_" + mod);
                        node.setProperty(PropertyType.LOCATION.toString(), matcher.group(0));
                        node.setProperty(PropertyType.TYPE.toString(), "p_Y");
                        node.setProperty(PropertyType.DISPLAY_NAME.toString(),  "p-Y_"+ matcher.group(0));
                        node.createRelationshipTo(current, RelTypes.PHOSPHORYLATION);
                    } else if (mod.contains("p-STY")) {
                        Node node = graphDb.createNode(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
                        node.setProperty(PropertyType.DB_ID.toString(), entityToString(entity) + "_" + mod);
                        node.setProperty(PropertyType.LOCATION.toString(), matcher.group(0));
                        node.setProperty(PropertyType.TYPE.toString(), "p_STY");
                        node.setProperty(PropertyType.DISPLAY_NAME.toString(),  "p-STY_"+ matcher.group(0));
                        node.createRelationshipTo(current, RelTypes.PHOSPHORYLATION);
                    }else{
                        Node node = graphDb.createNode(Label.label(LabelTypes.MODIFICATION.toString()));
                        Pattern p = Pattern.compile("\\[.*\\]");
                        Matcher m = p.matcher(mod);
                        m.find();
                        node.setProperty(PropertyType.DB_ID.toString(), entityToString(entity) + "_" + m.group(0));
                        node.setProperty(PropertyType.LOCATION.toString(), matcher.group(0));
                        node.setProperty(PropertyType.TYPE.toString(), m.group(0));
                        node.setProperty(PropertyType.DISPLAY_NAME.toString(),  m.group()+"_"+ matcher.group(0));
                        node.createRelationshipTo(current, RelTypes.MODIFICATION);
                    }
                }
            }
        } else {
            // it is a protein-PE
            // get its components uids
            Set<String> uids = getComponentStrings((Protein) entity);
            // if a there are UIDs
            creatUIDNode(uids);
        }
        return current;
    }


    private  void creatUIDNode(Set<String> uids){
        Pattern p = Pattern.compile(UID_PATTERN);

        if (!(uids.isEmpty())) {
            // iterate through (should usually only be one but comes in a set)
            for (String uid : uids) {
                // if updated
                if(uid.contains("_u")){
                    // get the base uid
                    Matcher m = p.matcher(uid);
                    if (m.find() ) {
                        String theGroup = m.group(0);
                        // if the UID node doesnt exist
                        if (graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), theGroup) == null) {
                            // make it and attach it to the Protein node
                            Node node = graphDb.createNode(Label.label(LabelTypes.UNIPROT_ID.toString()));
                            node.setProperty(PropertyType.UNIPROT_ID.toString(), theGroup);
                            node.setProperty(PropertyType.DB_ID.toString(), theGroup);
                            node.setProperty(PropertyType.DISPLAY_NAME.toString(), theGroup);
                            node.setProperty(PropertyType.TYPE.toString(), "ProteinID");
                            node.setProperty(PropertyType.STATUS.toString(), "Updated");
                            node.setProperty(PropertyType.DB_CONNECTION.toString(), ("https://www.uniprot.org/uniprot/" + theGroup));
                            node.createRelationshipTo(current, RelTypes.ID_BELONGS_TO);

                        } else {
                            graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), theGroup).createRelationshipTo(current, RelTypes.ID_BELONGS_TO);
                        }
                    }
                }else if(uid.contains("_d")){
                    // get the base uid
                    Matcher m = p.matcher(uid);
                    if (m.find() ) {
                        String theGroup = m.group(0);
                        // if the UID node doesnt exist
                        if (graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), theGroup) == null) {
                            // make it and attach it to the Protein node
                            Node node = graphDb.createNode(Label.label(LabelTypes.UNIPROT_ID.toString()));
                            node.setProperty(PropertyType.UNIPROT_ID.toString(), theGroup);
                            node.setProperty(PropertyType.DB_ID.toString(), theGroup);
                            node.setProperty(PropertyType.DISPLAY_NAME.toString(), theGroup);
                            node.setProperty(PropertyType.TYPE.toString(), "ProteinID");
                            node.setProperty(PropertyType.STATUS.toString(), "Deleted?");
                            node.createRelationshipTo(current, RelTypes.ID_BELONGS_TO);

                        } else {
                            graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), theGroup).createRelationshipTo(current, RelTypes.ID_BELONGS_TO);
                        }
                    }
                }else if(uid.contains("_c")){
                    // if the UID node doesnt exist
                    if (graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), uid) == null) {
                        // make it and attach it to the Protein node
                        Node node = graphDb.createNode(Label.label(LabelTypes.UNIPROT_ID.toString()));
                        node.setProperty(PropertyType.UNIPROT_ID.toString(), uid);
                        node.setProperty(PropertyType.DB_ID.toString(), uid);
                        node.setProperty(PropertyType.DISPLAY_NAME.toString(), uid);
                        node.setProperty(PropertyType.TYPE.toString(), "ProteinID");
                        node.setProperty(PropertyType.STATUS.toString(), "Current");
                        node.setProperty(PropertyType.DB_CONNECTION.toString(), ("https://www.uniprot.org/uniprot/" + uid));
                        node.createRelationshipTo(current, RelTypes.ID_BELONGS_TO);


                    } else {
                        graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), uid).createRelationshipTo(current, RelTypes.ID_BELONGS_TO);
                    }
                }else{
                    if (graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), uid) == null) {
                        // make it and attach it to the Protein node
                        Node node = graphDb.createNode(Label.label(LabelTypes.UNIPROT_ID.toString()));
                        node.setProperty(PropertyType.UNIPROT_ID.toString(), uid);
                        node.setProperty(PropertyType.DB_ID.toString(), uid);
                        node.setProperty(PropertyType.DISPLAY_NAME.toString(), uid);
                        node.setProperty(PropertyType.TYPE.toString(), "ProteinID");
                        node.setProperty(PropertyType.DB_CONNECTION.toString(), ("https://www.uniprot.org/uniprot/" + uid));
                        node.createRelationshipTo(current, RelTypes.ID_BELONGS_TO);


                    } else {
                        graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), uid).createRelationshipTo(current, RelTypes.ID_BELONGS_TO);
                    }

                }
            }
        }
    }

    /**
     * creates a physical entity node in the database
     */
    private   Node createPENode(Entity entity){

        // else if its  a PE thats not a prot just make the node INCLUDING COMPLEX
        current = graphDb.createNode(Label.label(classToString(entity.getClass().toString())));
        current.addLabel(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
        // UniqueID
        //take everything after the hashtag
        current.setProperty(PropertyType.DB_ID.toString(), entityToString(entity));
        if (entity.getDisplayName() == null) {
            current.setProperty(PropertyType.DISPLAY_NAME.toString(), "NA");
        } else {
            current.setProperty(PropertyType.DISPLAY_NAME.toString(), entity.getDisplayName());
        }
        String typeString = classToString(entity.getClass().toString());
        current.setProperty(PropertyType.TYPE.toString(), typeString);
        if (((PhysicalEntity) entity).getCellularLocation() == null) {
            current.setProperty(PropertyType.LOCATION.toString(), "NA");
        } else {
            current.setProperty(PropertyType.LOCATION.toString(), entityCellularLocation(entity));
        }
        Set<String> comments = entity.getComment();
        Set<String> uids = getComponentStrings((PhysicalEntity) entity);
        if (!(uids.isEmpty())) {
            HashSet<String> newUIDs = new HashSet<>();
            for (String uid : uids) {
                String newUID = "";
                if (uid.contains("#")) {
                    Pattern entPattern = Pattern.compile("(?<=#).*");
                    Matcher entMatcher = entPattern.matcher(uid);
                    entMatcher.find();
                    newUID = entMatcher.group(0);
                    newUIDs.add(newUID);
                }
                else{
                    newUIDs.add(uid);
                }
            }
            for(String comment: comments) {
                if(comment.contains("Converted from EntitySet")){
                    String uidList = String.join(", ", newUIDs);
                    current.setProperty(PropertyType.SET.toString(), uidList);
                }else{
                    String uidList = String.join(", ", newUIDs);
                    current.setProperty(PropertyType.UNIPROT_ID.toString(), uidList);
                }
            }
        }


        String dbconnection = databaseLink(entity);
        current.setProperty(PropertyType.DB_CONNECTION.toString(), dbconnection);

        return current;
    }

    /**
     * retrieves the UID's from the biopax file
     * @param entity
     * @return
     */
    private  Set<String> getUniprotID(Protein entity){
        // regex for a uniprot ID
        Pattern p = Pattern.compile(UID_PATTERN);
        Set<String> UID = new HashSet<String>();

        // if protein has an entityRef then get the xref from it
        EntityReference entityReference = entity.getEntityReference();
        // if there is an xref, get all
        if(!(entityReference == null)) {
            Set<Xref> xrefSet = entity.getEntityReference().getXref();
            for (Xref xref : xrefSet) {
                if(xref.toString().contains("UniProt") | xref.toString().contains("uniprot")) {
                    // if the xref is a uniprotID, get just the ID
                    Matcher m = p.matcher(xref.toString());
                    if (m.find() ) {
                        String theGroup = m.group(0);
                        // if uid and updated
                        if(xref.getDb().equals("UniProt_updated")){
                            UID.add(theGroup + "_u");
                            // if uid and deleted
                        }else if(xref.getDb().equals("UniProt_deleted")){
                            UID.add(theGroup + "_d");
                        }else if(xref.getDb().equals("UniProt_current")){
                            UID.add(theGroup + "_c");
                        }else{
                            UID.add(theGroup);
                        }
                    } else {
                        UID.add("No_UID_found");
                    }
                }
            }
        }
        return UID;
    }

    private  Set<PhysicalEntity> getComponents(PhysicalEntity entity, Set<PhysicalEntity> components){
        // Components -> if the components are complexes or have members themselves then recurse
        if (entity instanceof Complex){
            //System.out.println("complex: " + entity + " components " + ((Complex) entity).getComponent());
            Set<PhysicalEntity> entityComponentSet = ((Complex) entity).getComponent();
            //components.addAll(entityComponentSet);
            for(PhysicalEntity physicalEntity: entityComponentSet){
                if(physicalEntity instanceof Complex){
                    getComponents(physicalEntity,components);
                }else if(!(physicalEntity.getMemberPhysicalEntity().isEmpty())){
                    getComponents(physicalEntity, components);
                }else{
                    components.add(physicalEntity);
                }
            }
        }
        // Members -> if the members are complexes or have members themselves then recurse
        if(!(entity.getMemberPhysicalEntity().isEmpty())){
            //System.out.println("complex: " + entity + " members " + entity.getMemberPhysicalEntity());
            Set<PhysicalEntity> entityMemberSet = entity.getMemberPhysicalEntity();
            for(PhysicalEntity physicalEntity: entityMemberSet){
                if(physicalEntity instanceof Complex){
                    getComponents(physicalEntity,components);
                }else if(!(physicalEntity.getMemberPhysicalEntity().isEmpty())){
                    getComponents(physicalEntity, components);
                }else{
                    components.add(physicalEntity);
                }
            }
        }
        return components;
    }

    /**
     * code that sets up and calls getComponents and getUniprotID for things with multi UIDs
     * @param physicalEntity
     * @return
     */
    private  Set<String > getComponentStrings(PhysicalEntity physicalEntity){
        // empty component set
        Set<PhysicalEntity> physicalEntities = new HashSet<PhysicalEntity>();
        // full component set for complex
        Set<PhysicalEntity> components =  getComponents(physicalEntity, physicalEntities); // recursive
        // string list of components to add to
        Set<String> componentString = new HashSet<>();

        // if a component is a prot, get uid, else add string to list
        for(PhysicalEntity component: components){
            if (component instanceof Protein){
                Set<String> uids = getUniprotID((Protein) component);
                componentString.addAll(uids);
            }else{
                componentString.add(component.toString());
            }
        }
        return componentString;
    }


    private static ArrayList<String>  getMod(PhysicalEntity physicalEntity){
        ArrayList<String> mods = new ArrayList<>();

        // get all features for a pe
        Set<EntityFeature> features = physicalEntity.getFeature();
        // if not empty
        if(!features.isEmpty()) {
            for (EntityFeature feature : features) {
                if (feature.toString().contains("ModificationFeature: [O-phospho-L-serine]")) {
                    SequenceLocation featureLocation = feature.getFeatureLocation();
                    if (featureLocation == null){
                        mods.add("p-S_unknown");
                    }else {
                        mods.add("p-S_" +  feature.getFeatureLocation().toString());
                    }

                } else if (feature.toString().contains("ModificationFeature: [O-phospho-L-threonine]")) {
                    SequenceLocation featureLocation = feature.getFeatureLocation();
                    if (featureLocation == null){
                        mods.add("p-T_unknown");
                    }else {
                        mods.add("p-T_"+ feature.getFeatureLocation().toString());
                    }
                } else if (feature.toString().contains("ModificationFeature: [O4'-phospho-L-tyrosine]")) {
                    SequenceLocation featureLocation = feature.getFeatureLocation();
                    if (featureLocation == null){
                        mods.add("p-Y_unknown");
                    }else {
                        mods.add("p-Y_"+ feature.getFeatureLocation().toString());
                    }

                } else if (feature.toString().contains("ModificationFeature: [phosphorylated residue]")) {
                    SequenceLocation featureLocation = feature.getFeatureLocation();
                    if (featureLocation == null){
                        mods.add("p-STY_unknown");
                    }else {
                        mods.add("p-STY_"+ feature.getFeatureLocation().toString());
                    }

                }else {
                    if (feature.toString().contains("ModificationFeature: [")) {
                        SequenceLocation featureLocation = feature.getFeatureLocation();
                        if (featureLocation == null) {
                            mods.add(feature.toString() + "_unknown");
                        } else {
                            mods.add(feature.toString() + "_" + featureLocation.toString());
                        }
                    }
                }

            }
            // if features is empty, check for mods on its ent ref
        }else{
            //HPRD SPECIFIC - mods are on entity refs
            EntityReference entityReference = ((Protein) physicalEntity).getEntityReference();
            Set<EntityFeature> entityFeatures = entityReference.getEntityFeature();
            if(!entityFeatures.isEmpty()){
                for (EntityFeature feature : entityFeatures) {
                    if (feature.toString().contains("Phosphoserine")) {
                        SequenceLocation featureLocation = feature.getFeatureLocation();
                        if (featureLocation == null){
                            mods.add("p-S_unknown");
                        }else {
                            mods.add("p-S_" +  feature.getFeatureLocation().toString());
                        }

                    } else if (feature.toString().contains("Phosphothreonine")) {
                        SequenceLocation featureLocation = feature.getFeatureLocation();
                        if (featureLocation == null){
                            mods.add("p-T_unknown");
                        }else {
                            mods.add("p-T_"+ feature.getFeatureLocation().toString());
                        }
                    } else if (feature.toString().contains("Phosphotyrosine")) {
                        SequenceLocation featureLocation = feature.getFeatureLocation();
                        if (featureLocation == null){
                            mods.add("p-Y_unknown");
                        }else {
                            mods.add("p-Y_"+ feature.getFeatureLocation().toString());
                        }

                    }else {
                        if (feature.toString().contains("ModificationFeature: [")) {
                            SequenceLocation featureLocation = feature.getFeatureLocation();
                            if (featureLocation == null) {
                                mods.add(feature.toString() + "_unknown");
                            } else {
                                mods.add(feature.toString() + "_" + featureLocation.toString());
                            }
                        }
                    }

                }
            }
        }
        return mods;
    }

    /**
     * converts the entity to a string
     * @param entity
     * @return
     */
    private String entityToString( Entity entity){
        String entString = entity.toString();

        if (entity.toString().contains(xmlBase)){
            String entityString = entity.toString();
            Pattern entPattern = Pattern.compile("(?<=#).*");
            Matcher entMatcher = entPattern.matcher(entityString);
            entMatcher.find();
            entString = entMatcher.group(0);
        }

        return entString;
    }


    /**
     * e.g.
     * takes a long string 'class org.iiusdfhb.level3.ProteinImpl' and turns it into just Protein
     * @param classString
     * @return
     */
    private static String classToString(String classString){

        Matcher m = Pattern.compile("3\\.(.*)Impl$").matcher(classString);
        if (m.find()) {
            classString = m.group(1);
        }
        return classString;
    }

    /**
     * Takes an entity and find the stable reactome identfier (website Link) and links it
     * @param entity
     * @return
     */
    private String databaseLink(Entity entity) {
        String dbconnecton = "";
        for (Xref xref : entity.getXref()) {
            for (String s : xref.getComment()) {
                if(s.startsWith("Reactome stable identifier.")){
                    Matcher m = Pattern.compile(":(.*)").matcher(s);
                    if (m.find()) {
                        dbconnecton = m.group(1);
                    }
                }
            }
        }
        return dbconnecton;
    }

    /**
     * takes the location string from model
     * @param entity
     * @return
     */
    private static String entityCellularLocation(Entity entity){
        String peCL = ((PhysicalEntity) entity).getCellularLocation().toString();
        int i = peCL.lastIndexOf('_');
        peCL = peCL.substring(i+1);

        return peCL;
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

    /**
     * Will take a biopax model and tag all UniProt xrefs with 'updated' or 'deleted'
     * It will not overwrite the file
     * @param model
     * the given BioPax Model. Assumes UniProt DB xrefs are "UniProt" or "uniprot"
     * @return a model where the UniProt DB xrefs are tagged to be read by createDBfromOWL
     *
     */
    Model updateBiopaxModel( Model model){
        // if not human may still want to keep (eg HIV)

        StopWatch stopwatch = new StopWatch();
        stopwatch.start();
        System.out.println("Updating UniProt ID's live ... \nCurrently downloading 2 UniProt files... ");


        // get secondary upniprot id's live and read into a dict
        HashMap<String, String> sec_ac = null;
        URL sec_acc_url = null;
        try {
            sec_acc_url = new URL(secondaryAccessionFile);
        } catch (MalformedURLException e) {
            e.printStackTrace();
            System.exit(1);
        }

        // all human uids currently in uniprot
        HashSet<String> allUIDs = new HashSet<>();
        URL IDs = null;
        String specStr = "";
        if(species.equalsIgnoreCase("h") | species.equalsIgnoreCase("human") | species.equals("9606")){
            specStr = "Human";
            try {
                IDs = new URL(humanUIDs);
            } catch (MalformedURLException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }else if(species.equalsIgnoreCase("m") | species.equalsIgnoreCase("mouse") | species.equals("10090")){
            specStr = "Mouse";
            try {
                IDs = new URL(mouseUIDs);
            } catch (MalformedURLException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }


        try {
            sec_ac = secondaryAccessionDictFromFile(sec_acc_url);
            System.out.println("Secondary Accession file downloaded");
            allUIDs = getCurrentUIDs(IDs);
            System.out.println("Current UniProt ID file downloaded");
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }


        Integer updatedCounter = 0;
        Integer delOrNotHumanCounter = 0;


        Set<Xref> xrefs = model.getObjects(Xref.class);
        Pattern p = Pattern.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}");
        Pattern p2 = Pattern.compile("(\\-[0-9*]{1,2})");


        for(Xref xref: xrefs){
            if(xref.toString().contains("UniProt") | xref.toString().contains("uniprot")) {
                Matcher m = p.matcher(xref.toString());
                m.find();
                String UID = m.group(0);

                // if it's an isoform - go through normal and update/delete but add iso number to it
                Matcher m2 = p2.matcher(xref.toString());
                if(m2.find()){
                    String isoNum = m2.group(0);
                    // if its a secondary accession update its xref
                    if (sec_ac.containsKey(UID)) {
                        String newUID = sec_ac.get(UID);
                        // check to make sure the new uid is human
                        if (allUIDs.contains(newUID)) {
                            xref.setId(newUID + isoNum);
                            xref.setDb("UniProt_updated" );
                            updatedCounter++;
                        }else{
                            xref.setDb("UniProt_deleted");
                            delOrNotHumanCounter ++;
                        }
                        // else its primary check if its not current human and tag it
                    } else if (!allUIDs.contains(UID)) {
                        // its deleted
                        xref.setDb("UniProt_deleted");
                        delOrNotHumanCounter ++;
                    }else{
                        xref.setDb("UniProt_current");
                    }

                    // if its not an isoform
                    // if its a secondary accession update its xref
                }else if (sec_ac.containsKey(UID)) {
                    String newUID = sec_ac.get(UID);
                    // check to make sure the new uid is human
                    if (allUIDs.contains(newUID)) {
                        xref.setId(newUID);
                        xref.setDb("UniProt_updated" );
                        updatedCounter++;
                    }else{
                        xref.setDb("UniProt_deleted");
                        delOrNotHumanCounter ++;
                    }
                    // else its not an isoform
                    // it is a primary acc check if its not current human and tag it
                } else if (!allUIDs.contains(UID)) {
                    // its deleted
                    xref.setDb("UniProt_deleted");
                    delOrNotHumanCounter ++;
                }else{
                    xref.setDb("UniProt_current");
                }
            }
        }

        System.out.println("Number of ID's updated:"+ updatedCounter);
        System.out.println("Number of ID's flagged as deleted or not "+ specStr+": " + delOrNotHumanCounter);
        System.out.println(stopwatch.toString());
        stopwatch.stop();

        return model;
    }


    /**
     * takes a url with a file of secondary accs mapped to primary accs and returns a hashmap of the mapping
     * @param url
     * @return HashMap with secondary accessions as keys and primary accessions as values
     * @throws IOException for input url
     */
    HashMap<String, String> secondaryAccessionDictFromFile(URL url ) throws IOException{
        HashMap<String, String> secAccessions = new HashMap<>();

        BufferedReader BR = new BufferedReader(new InputStreamReader(url.openStream()));
        String line = "";

        Pattern p = Pattern.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}");

        while ((line = BR.readLine()) != null) {

            Matcher m = p.matcher(line);
            String secondary = "";
            String primary = "";

            if (m.find()) {
                String[] splitLine = line.split("\\s\\s\\s\\s*");
                secondary = splitLine[0];
                primary = splitLine[1];

                secAccessions.put(secondary, primary);

            }
        }
        return secAccessions;
    }

    /**
     *
     * @param url
     * should be the current version of UniProt
     * @return HashSet of current UniProt ID's found in given url. One per line.
     * @throws IOException
     * for input url
     */
    HashSet<String> getCurrentUIDs(URL url) throws IOException {
        //TODO if file already exists use that if it's date of creation is over 1 month old or it doesnt exist download it
        BufferedReader BR = new BufferedReader(new InputStreamReader(url.openStream()));

        HashSet<String> currentUIDset = new HashSet<>();
        String line = BR.readLine();

        while((line = BR.readLine())!= null){
            currentUIDset.add(line);
        }
        return currentUIDset;
    }

    HashSet<String> getLabelledUIDs(String urlString) throws IOException {

        URL url = null;

        try {
            url = new URL(urlString);
        } catch (MalformedURLException e) {
            e.printStackTrace();
        }

        BufferedReader BR = new BufferedReader(new InputStreamReader(url.openStream()));

        HashSet<String> set = new HashSet<>();
        String line ="";

        while((line = BR.readLine())!= null){
            set.add(line);

        }

        return set;
    }

    void addLabelsToDB(String property, String urlString){

        HashSet<String> labelledUIDs = new HashSet<>();
        try {
            labelledUIDs = getLabelledUIDs(urlString);
        } catch (IOException e) {
            System.err.println("Could not read in file from URL: " + urlString);
            e.printStackTrace();
            System.exit(1);
        }

        try (Transaction tx = graphDb.beginTx()) {

            for (String uid: labelledUIDs) {
                Node node = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), uid);

                if(node != null){
                    Iterable<Relationship> relationships = node.getRelationships(RelTypes.ID_BELONGS_TO);
                    for(Relationship rel: relationships){
                        Node endNode = rel.getEndNode();
                        endNode.setProperty(property, property);
                    }
                }
            }
            tx.success();
        }
    }

    public void addUniProtGeneNames() {

        URL IDs = null;
        String specStr = "";
        if(species.equalsIgnoreCase("h") | species.equalsIgnoreCase("human") | species.equals("9606")){
            specStr = "Human";
            try {
                IDs = new URL(humanUniProtNames);
            } catch (MalformedURLException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }else if(species.equalsIgnoreCase("m") | species.equalsIgnoreCase("mouse") | species.equals("10090")){
            specStr = "Mouse";
            try {
                IDs = new URL(mouseUniProtNames);
            } catch (MalformedURLException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        try(Transaction tx = graphDb.beginTx()){

            try{
                BufferedReader BR = new BufferedReader(new InputStreamReader(IDs.openStream()));
                String line = BR.readLine();
                while((line = BR.readLine())!= null){
                    String[] split = line.split("\t");
                    if(split.length == 2){
                        String UID = split[0];
                        String name = split[1];

                        Node uidNode = graphDb.findNode(Label.label(LabelTypes.UNIPROT_ID.toString()), PropertyType.UNIPROT_ID.toString(), UID);
                        if(uidNode != null){
                            if(!uidNode.hasProperty(PropertyType.UNIPROT_NAME.toString())){
                                uidNode.setProperty(PropertyType.UNIPROT_NAME.toString(), name);

                                Iterable<Relationship> relationships = uidNode.getRelationships(RelTypes.ID_BELONGS_TO);
                                for(Relationship rel: relationships){
                                    Node endNode = rel.getEndNode();
                                    if(endNode.hasProperty(PropertyType.UNIPROT_NAME.toString())){
                                        String s = endNode.getProperty(PropertyType.UNIPROT_NAME.toString()).toString();
                                        s = s.concat(", "+ name);
                                        endNode.setProperty(PropertyType.UNIPROT_NAME.toString(), s);
                                    }else{
                                        endNode.setProperty(PropertyType.UNIPROT_NAME.toString(),name);
                                    }

                                }
                            }
                        }

                    }
                }
            }catch (IOException e){
                System.out.println("Could not download file: " + IDs);
                e.printStackTrace();
                System.exit(1);
            }
            tx.success();
        }


    }

    Boolean relationshipBetween(Node n1, Node n2) { // RelationshipType type, Direction direction
        for (Relationship rel : n1.getRelationships()) { // n1.getRelationships(type,direction)
            if (rel.getOtherNode(n1).equals(n2)){
                return true;
            }
        }
        return false;
    }

}

