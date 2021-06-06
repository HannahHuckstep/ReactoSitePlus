package org.wehi.hucksteph;

import org.apache.commons.io.FileUtils;
import org.biopax.paxtools.model.BioPAXFactory;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Protein;
import org.biopax.paxtools.model.level3.UnificationXref;
import org.biopax.paxtools.model.level3.Xref;
import org.junit.jupiter.api.Test;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.unsafe.impl.batchimport.input.InputException;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class DatabaseFactoryTest  {

    private final File TEST_OWL_FILE = new File("test/InsulinTestNtwk.owl");
    private final File OTHER_TEST_OWL_FILE = new File("test/OtherTestNtwk.owl");
    private final File DATABASE_ACTUAL_PATH = new File("test/actual/");

    @Test
    public void dbNotEmptyTest() {

        File tempDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraph = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraph, false, "h");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original
        // just made
        GraphDatabaseService actual = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraph);

        // cypher query that returns everything
        String a = null;
        try(Transaction tx = actual.beginTx()){
            Result result = actual.execute("MATCH (n) RETURN n");
            a = result.resultAsString();
            tx.success();
        }

        String e = null;


        assertNotEquals(e,a);

        try{
            FileUtils.deleteDirectory(tempDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    @Test
    public void dbShouldNotUpdateTest(){
        File tempDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraph = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");

        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraph, false, "h");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original
        // just made
        GraphDatabaseService actual = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraph);

        try(Transaction tx = actual.beginTx()){
            ResourceIterator<Node> nodes = actual.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            for (ResourceIterator<Node> it = nodes; it.hasNext(); ) {
                Node node = it.next();
                if(node.hasProperty(PropertyType.STATUS.toString())){
                    fail();
                }
            }
            tx.success();
        }

        try{
            FileUtils.deleteDirectory(tempDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
    }

    @Test
    public void dbShouldUpdateDefaultTest(){
        File tempDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraph = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraph);
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original
        // just made
        GraphDatabaseService actual = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraph);

        try(Transaction tx = actual.beginTx()){
            ResourceIterator<Node> nodes = actual.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            while (nodes.hasNext()){
                Node next = nodes.next();
                if (next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("Q9Y6L4")){
                    assertEquals("Updated", next.getProperty(PropertyType.STATUS.toString()).toString());
                }else if(next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("U3BNT9")){
                    assertEquals("Deleted?", next.getProperty(PropertyType.STATUS.toString()).toString());
                }else if(next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("P06213")){
                    assertEquals("Current", next.getProperty(PropertyType.STATUS.toString()).toString());
                }else if(next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("A5A626")){
                    assertEquals("Deleted?", next.getProperty(PropertyType.STATUS.toString()).toString());
                }else if(next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("Q9UH92-2")){
                    assertEquals("Current", next.getProperty(PropertyType.STATUS.toString()).toString());
                }
            }

            tx.success();
        }
        try{
            FileUtils.deleteDirectory(tempDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
    }

    @Test
    public void dbShouldUpdateTest(){
        File tempDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraph = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraph, true, "h");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original
        // just made
        GraphDatabaseService actual = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraph);

        try(Transaction tx = actual.beginTx()){
            ResourceIterator<Node> nodes = actual.findNodes(Label.label(LabelTypes.UNIPROT_ID.toString()));
            while (nodes.hasNext()){
                Node next = nodes.next();
                if (next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("Q9Y6L4")){
                    assertEquals("Updated", next.getProperty(PropertyType.STATUS.toString()).toString());
                }else if(next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("U3BNT9")){
                    assertEquals("Deleted?", next.getProperty(PropertyType.STATUS.toString()).toString());
                }else if(next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("P06213")){
                    assertEquals("Current", next.getProperty(PropertyType.STATUS.toString()).toString());
                }else if(next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("A5A626")){
                    assertEquals("Deleted?", next.getProperty(PropertyType.STATUS.toString()).toString());
                }else if(next.getProperty(PropertyType.UNIPROT_ID.toString()).equals("Q9UH92-2")){
                    assertEquals("Current", next.getProperty(PropertyType.STATUS.toString()).toString());
                }
            }

            tx.success();
        }
        try{
            FileUtils.deleteDirectory(tempDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    @Test
    public void updateDeletedUIDTest(){
        BioPAXFactory factory = BioPAXLevel.L3.getDefaultFactory();
        Model model = factory.createModel();
        Protein P1 = model.addNew(Protein.class, "Protein1");
        UnificationXref uxref = model.addNew(UnificationXref.class,"huckstep.wehi.org#xref_A5A626");
        uxref.setDb("UniProt");
        uxref.setId("A5A626");
        P1.addXref(uxref);

        DatabaseFactory dbf = new DatabaseFactory("h");
        Model updatedBiopaxModel = dbf.updateBiopaxModel(model);

        Set<Xref> uodatedXrefSet = updatedBiopaxModel.getObjects(Xref.class);
        String actual = "";
        for(Xref updatedXref: uodatedXrefSet){
            actual = updatedXref.toString();
        }

        assertEquals("UniProt_deleted:A5A626", actual);
    }

    @Test
    public void updateSecondaryUIDTest(){
        BioPAXFactory factory = BioPAXLevel.L3.getDefaultFactory();
        Model model = factory.createModel();
        Protein P1 = model.addNew(Protein.class, "Protein1");
        UnificationXref uxref = model.addNew(UnificationXref.class,"huckstep.wehi.org#xref_A5A626");
        uxref.setDb("UniProt");
        uxref.setId("Q9Y6L4");
        P1.addXref(uxref);

        DatabaseFactory dbf = new DatabaseFactory("h");
        Model updatedBiopaxModel = dbf.updateBiopaxModel(model);

        Set<Xref> uodatedXrefSet = updatedBiopaxModel.getObjects(Xref.class);
        String actual = "";
        for(Xref updatedXref: uodatedXrefSet){
            actual = updatedXref.toString();
        }

        assertEquals("UniProt_updated:Q92838", actual);
    }

    @Test
    public void updateNonHumanSecondaryUIDTest(){
        BioPAXFactory factory = BioPAXLevel.L3.getDefaultFactory();
        Model model = factory.createModel();
        Protein P1 = model.addNew(Protein.class, "Protein1");
        UnificationXref uxref = model.addNew(UnificationXref.class,"huckstep.wehi.org#xref_A5A626");
        uxref.setDb("UniProt");
        uxref.setId("U3BNT9");
        P1.addXref(uxref);

        DatabaseFactory dbf = new DatabaseFactory("h");
        Model updatedBiopaxModel = dbf.updateBiopaxModel(model);

        Set<Xref> uodatedXrefSet = updatedBiopaxModel.getObjects(Xref.class);
        String actual = "";
        for(Xref updatedXref: uodatedXrefSet){
            actual = updatedXref.toString();
        }

        assertEquals("UniProt_deleted:U3BNT9", actual);

    }

    @Test
    public void dontUpdateCurrentUIDTest(){
        BioPAXFactory factory = BioPAXLevel.L3.getDefaultFactory();
        Model model = factory.createModel();
        Protein P1 = model.addNew(Protein.class, "Protein1");
        UnificationXref uxref = model.addNew(UnificationXref.class,"huckstep.wehi.org#xref_A5A626");
        uxref.setDb("UniProt");
        uxref.setId("P06213");
        P1.addXref(uxref);

        DatabaseFactory dbf = new DatabaseFactory("h");
        Model updatedBiopaxModel = dbf.updateBiopaxModel(model);

        Set<Xref> uodatedXrefSet = updatedBiopaxModel.getObjects(Xref.class);
        String actual = "";
        for(Xref updatedXref: uodatedXrefSet){
            actual = updatedXref.toString();
        }

        assertEquals("UniProt_current:P06213", actual);

    }

    @Test
    void updatingIsoformTest(){
        BioPAXFactory factory = BioPAXLevel.L3.getDefaultFactory();
        Model model = factory.createModel();
        Protein P1 = model.addNew(Protein.class, "Protein1");
        UnificationXref uxref = model.addNew(UnificationXref.class,"huckstep.wehi.org#xref_B2RAV8-2");
        uxref.setDb("UniProt");
        uxref.setId("B2RAV8-2");
        P1.addXref(uxref);

        DatabaseFactory dbf = new DatabaseFactory("h");
        Model updatedBiopaxModel = dbf.updateBiopaxModel(model);

        Set<Xref> uodatedXrefSet = updatedBiopaxModel.getObjects(Xref.class);
        String actual = "";
        for(Xref updatedXref: uodatedXrefSet){
            actual = updatedXref.toString();
        }

        assertEquals("UniProt_updated:Q9UH92-2", actual);

    }

    @Test
    void updateMouseUIDTest(){
        BioPAXFactory factory = BioPAXLevel.L3.getDefaultFactory();
        Model model = factory.createModel();
        Protein P1 = model.addNew(Protein.class, "Protein1");
        UnificationXref uxref = model.addNew(UnificationXref.class,"huckstep.wehi.org#xref_B2RAV8-2");
        uxref.setDb("UniProt");
        uxref.setId("P15208");
        P1.addXref(uxref);

        DatabaseFactory dbf = new DatabaseFactory("m");
        Model updatedBiopaxModel = dbf.updateBiopaxModel(model);

        Set<Xref> uodatedXrefSet = updatedBiopaxModel.getObjects(Xref.class);
        String actual = "";
        for(Xref updatedXref: uodatedXrefSet){
            actual = updatedXref.toString();
        }

        assertEquals("UniProt_current:P15208", actual);

    }

    @Test
    void getKinases() {
        String KINASES = "https://www.uniprot.org/uniprot/?query=keyword:%22Kinase%20[KW-0418]%22&format=list&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22";
        DatabaseFactory dbf = new DatabaseFactory();
        HashSet<String> hashSet = new HashSet<>();
        try {
            hashSet = dbf.getLabelledUIDs(KINASES);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // testing to see it contains randomly selected kinases as list may be added to or removed from in future
        if(!hashSet.contains("Q9Y3S1")){
            fail();
        }else if (!hashSet.contains("O60285")){
            fail();
        }else if (!hashSet.contains("Q9NQU5")){
            fail();
        }else if (!hashSet.contains("Q86UX6")){
            fail();
        }else if (!hashSet.contains("Q86TB3")){
            fail();
        }else if (!hashSet.contains("A2RU49")){
            fail();
        }

    }

    @Test
    void getTranscriptionFactors() {
        String TRANSCRIPTION_FACTORS = "https://www.uniprot.org/uniprot/?query=goa:(%22DNA-binding%20transcription%20factor%20activity%20[3700]%22)%20(reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22)&format=list";
        DatabaseFactory dbf = new DatabaseFactory();
        HashSet<String> hashSet = new HashSet<>();
        try {
            hashSet = dbf.getLabelledUIDs(TRANSCRIPTION_FACTORS);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // testing to see it contains randomly selected kinases as list may be added to or removed from in future
        if(!hashSet.contains("Q9UMR3")){
            fail();
        }else if (!hashSet.contains("Q52M93")){
            fail();
        }else if (!hashSet.contains("O14813")){
            fail();
        }else if (!hashSet.contains("Q5JVG8")){
            fail();
        }else if (!hashSet.contains("P15923")){
            fail();
        }else if (!hashSet.contains("P21506")){
            fail();
        }
    }

    @Test
    void getCellSurfaceReceptors() {
        String CELL_SURFACE_RECEPTORS = "https://www.uniprot.org/uniprot/?query=goa:(%22cell%20surface%20receptor%20signaling%20pathway%20involved%20in%20cell-cell%20signaling%20[1905114]%22)%20(reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22)&format=list";
        DatabaseFactory dbf = new DatabaseFactory();
        HashSet<String> hashSet = new HashSet<>();
        try {
            hashSet = dbf.getLabelledUIDs(CELL_SURFACE_RECEPTORS);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // testing to see it contains randomly selected kinases as list may be added to or removed from in future
        if(!hashSet.contains("Q9NQB0")){
            fail();
        }else if (!hashSet.contains("P49840")){
            fail();
        }else if (!hashSet.contains("O00233")){
            fail();
        }else if (!hashSet.contains("P63208")){
            fail();
        }else if (!hashSet.contains("Q9BXB1")){
            fail();
        }else if (!hashSet.contains("O15130")){
            fail();
        }

    }

    @Test
    void extraLabelsTest(){

        File tempOutputDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraphDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");

        // in case graph wasn't already deleted
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        // Make a new graph
        DatabaseFactory dbf = new DatabaseFactory(OTHER_TEST_OWL_FILE, tempGraphDir, false, "human");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original

        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()) {
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));

            while (proteins.hasNext()){
                Node prot = proteins.next();

                if(prot.getProperty(PropertyType.DISPLAY_NAME.toString()).equals("Homologues of 2xHC-INS(25-54)")){ // Q9Y3S1 KINASE, Q9NQB0 TF & Cell Surface
                    if(!prot.hasProperty(PropertyType.KINASE.toString())){
                        fail();
                    }else if(!prot.hasProperty(PropertyType.TRANSCRIPTION_FACTOR.toString())){
                        fail();
                    }else if(!prot.hasProperty(PropertyType.CELL_SURFACE_RECEPTOR.toString())){
                        fail();
                    }
                }else if(prot.getProperty(PropertyType.DISPLAY_NAME.toString()).equals("P25799")){
                    if(!prot.hasProperty(PropertyType.CELL_SURFACE_RECEPTOR.toString())){
                        fail();
                    }
                }
            }

            tx.success();
        }

        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
    }

    @Test
    void UniProtGeneNames() {

        File tempOutputDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraphDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");

        // in case graph wasn't already deleted
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        // Make a new graph
        DatabaseFactory dbf = new DatabaseFactory(OTHER_TEST_OWL_FILE, tempGraphDir, false, "human");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original


        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()) {
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));

            while (proteins.hasNext()){
                Node prot = proteins.next();

                if(prot.getProperty(PropertyType.DISPLAY_NAME.toString()).equals("Homologues of 4xHC-INS(90-110)")){ // Q9Y3S1 KINASE, Q9NQB0 TF & Cell Surface
                    String a = prot.getProperty(PropertyType.UNIPROT_NAME.toString()).toString();
                    assertEquals("WNK2, TCF7L2", a);
                }else if(prot.getProperty(PropertyType.DISPLAY_NAME.toString()).equals("P25799")){
                    String a = prot.getProperty(PropertyType.UNIPROT_NAME.toString()).toString();
                    assertEquals("TSPAN12", a);
                } else if(prot.getProperty(PropertyType.DISPLAY_NAME.toString()).equals("P06213")){
                    String a = prot.getProperty(PropertyType.UNIPROT_NAME.toString()).toString();
                    assertEquals("INSR", a);
                }
            }

            tx.success();
        }

        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
    }

    @Test
    void addingPathways(){

        File tempOutputDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraphDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");

        // in case graph wasn't already deleted
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        // Make a new graph
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraphDir, false, "human");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original


        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()) {
            ResourceIterator<Node> pathways = graphDb.findNodes(Label.label(LabelTypes.PATHWAY.toString()));
            while(pathways.hasNext()){
                Node pth = pathways.next();
                assertEquals("R-HSA-5660489", pth.getProperty(PropertyType.DB_ID.toString()));
                assertEquals("Hannah's Insulin Test Pathway", pth.getProperty(PropertyType.DISPLAY_NAME.toString()));
            }
            tx.success();
        }

        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
    }

    @Test
    public void throwErrorWhenNotHumanOrMouse(){
        File tempOutputDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraphDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");

        // in case graph wasn't already deleted
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        // Make a new graph
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraphDir, false, "THROW_ERROR");
        Exception exception = assertThrows(InputException.class, () ->
        {dbf.createDBfromOWL();
        });
        System.setIn(sysInBackup);// reset System.in to its original


        String expectedMessage = "Species must be either Human";
        String actualMessage = exception.getMessage();

        assertTrue(actualMessage.contains(expectedMessage));

        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

}