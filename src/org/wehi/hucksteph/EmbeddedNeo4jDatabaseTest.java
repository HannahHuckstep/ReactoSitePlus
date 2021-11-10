package org.wehi.hucksteph;

import org.apache.commons.io.FileUtils;
import org.biopax.paxtools.model.level2.complex;
import org.junit.jupiter.api.Test;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

import java.io.*;
import java.util.HashSet;

import static org.junit.jupiter.api.Assertions.*;

class EmbeddedNeo4jDatabaseTest {

    private final File TEST_OWL_FILE = new File("test/InsulinTestNtwk.owl");
    private final File TEST_DATA_FILE = new File("test/mapTestData.txt");
    private final File DATABASE_EXPECTED_PATH = new File("test/expected/");
    private final File DATABASE_ACTUAL_PATH = new File("test/actual/");

    @Test
    void testPrintDatabase() {
        File tempDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraph = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraph);
        System.setIn(sysInBackup);// reset System.in to its original

        try(Transaction tx = graphDb.beginTx()) {
            Node p1 = graphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            p1.setProperty(PropertyType.DISPLAY_NAME.toString(), "p1");

            Node p2 = graphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            p2.setProperty(PropertyType.DISPLAY_NAME.toString(), "p2");

            Node cplx1 = graphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            cplx1.setProperty(PropertyType.DISPLAY_NAME.toString(), "cplx1");

            Node ph1 = graphDb.createNode(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
            ph1.setProperty(PropertyType.DISPLAY_NAME.toString(), "ph1");

            Node uid1 = graphDb.createNode(Label.label(LabelTypes.UNIPROT_ID.toString()));
            uid1.setProperty(PropertyType.DISPLAY_NAME.toString(), "uid1");

            Node rxn1 = graphDb.createNode(Label.label(LabelTypes.INTERACTION.toString()));
            rxn1.setProperty(PropertyType.DISPLAY_NAME.toString(), "rxn1");

            p1.createRelationshipTo(rxn1, RelTypes.INPUT);
            p2.createRelationshipTo(rxn1, RelTypes.INPUT);
            cplx1.createRelationshipTo(rxn1, RelTypes.OUTPUT);
            p1.createRelationshipTo(cplx1, RelTypes.COMPONENT);
            p2.createRelationshipTo(cplx1, RelTypes.COMPONENT);
            ph1.createRelationshipTo(p1, RelTypes.PHOSPHORYLATION);
            uid1.createRelationshipTo(p1, RelTypes.ID_BELONGS_TO);

            tx.success();
        }
        graphDb.shutdown();

        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraph, tempDir);

        //capture std out
        // Create a stream to hold the output
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(baos);
        // IMPORTANT: Save the old System.out!
        PrintStream old = System.out;
        // Tell Java to use your special stream
        System.setOut(ps);
        // Print some output: goes to your special stream
        edb.printDatabase();
        // Put things back
        System.out.flush();
        System.setOut(old);
        String actual = baos.toString();


        try{
            FileUtils.deleteDirectory(tempDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

        String expected = "\n" +
                "Nodes:\n" +
                "\n" +
                "{DISPLAY_NAME=p1}\n" +
                "{DISPLAY_NAME=p2}\n" +
                "{DISPLAY_NAME=cplx1}\n" +
                "{DISPLAY_NAME=ph1}\n" +
                "{DISPLAY_NAME=uid1}\n" +
                "{DISPLAY_NAME=rxn1}\n" +
                "\n" +
                "Relationships:\n" +
                "\n" +
                "p1 -> INPUT{} -> rxn1\n" +
                "p2 -> INPUT{} -> rxn1\n" +
                "cplx1 -> OUTPUT{} -> rxn1\n" +
                "p1 -> COMPONENT{} -> cplx1\n" +
                "p2 -> COMPONENT{} -> cplx1\n" +
                "ph1 -> PHOSPHORYLATION{} -> p1\n" +
                "uid1 -> ID_BELONGS_TO{} -> p1\n";

        assertEquals(expected, actual);


        //System.out.println(outContent.toString());

    }

    @Test
    void testPrintLabelNumber() {

        File tempDir = new File(DATABASE_EXPECTED_PATH+ "/toBeDeleted/");
        File tempGraph = new File(DATABASE_EXPECTED_PATH+ "/toBeDeleted/GRAPH/");

        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraph);

        try(Transaction tx = graphDb.beginTx()){
            Node p1 = graphDb.createNode(Label.label("lb1"));
            Node p2 = graphDb.createNode(Label.label("lb2"));
            Node p3 = graphDb.createNode(Label.label("lb2"));
            Node p4 = graphDb.createNode(Label.label("lb3"));
            tx.success();
        }
        graphDb.shutdown();

        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraph, tempDir);

        //capture std out
        // Create a stream to hold the output
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(baos);
        // IMPORTANT: Save the old System.out!
        PrintStream old = System.out;
        // Tell Java to use your special stream
        System.setOut(ps);
        // Print some output: goes to your special stream
        edb.printLabelNumber("lb2");
        // Put things back
        System.out.flush();
        System.setOut(old);
        String actual = baos.toString();


        try{
            FileUtils.deleteDirectory(tempDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

        String expected = "There are 2 lb2's in this database: test/expected/toBeDeleted/GRAPH\n";

        assertEquals(expected, actual);

    }

    @Test
    void testMapMQPhosphopeps_HighestSupport_abundance() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "HighestSupport"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if proteins have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(12.9, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(5.95, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(3.2, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_0")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_0")){
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
    void testMapMQPhosphopeps_Median_abundance() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "Median"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(8.7, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(8.7, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(8.7, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_0")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_0")){
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
    void testMapMQPhosphopeps_Max_abundance() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "Max"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(16.5, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(16.5, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(16.5, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_0")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_0")){
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
    void testMapMQPhosphopeps_Mean_abundance() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "Mean"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(8.82, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(8.82, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(8.82, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_0")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_0")){
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
    void testMapMQPhosphopeps_Extreme_abundance() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "Extreme"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(16.5, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(16.5, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(16.5, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_0")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_0")){
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
    void testMapMQPhosphopeps_HighestSupport_lfc() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "HighestSupport"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(1.2, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(0.51, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(-0.08, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_1")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_1")){
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
    void testMapMQPhosphopeps_Median_lfc() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "Median"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(-0.082, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(-0.082, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(-0.082, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_1")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_1")){
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
    void testMapMQPhosphopeps_Max_lfc() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "Max"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(3.4, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(3.4, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(3.4, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_1")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_1")){
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
    void testMapMQPhosphopeps_Mean_lfc() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "Mean"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(0.584, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(0.584, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(0.584, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_1")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_1")){
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
    void testMapMQPhosphopeps_Extreme_lfc() {

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "Extreme"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> proteins = graphDb.findNodes(Label.label("Protein"));
            while(proteins.hasNext()){
                Node prot = proteins.next();
                if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("INSR(763-1382)")){ // unmodified
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.78, Double.valueOf(support_score_0));
                    assertEquals(3.4, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-3-INSR(400-1382)")){ // 3 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(1.08, Double.valueOf(support_score_0));
                    assertEquals(3.4, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("phospho-p-2-INSR(28-758)")){ // 2 phosns
                    String support_score_0 = prot.getProperty("SUPPORT_SCORE_1").toString();
                    String abundance_score_0 = prot.getProperty("ABUNDANCE_SCORE_1").toString();
                    assertEquals(0.98, Double.valueOf(support_score_0));
                    assertEquals(3.4, Double.valueOf(abundance_score_0));
                }else if((prot.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("4xHC-INS(90-110)")){ // not mapped
                    if(prot.hasProperty("SUPPORT_SCORE_1")){
                        fail();
                    }else if(prot.hasProperty("ABUNDANCE_SCORE_1")){
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
    void testScoreComplexes(){

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "HighestSupport"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }

        // text to see if complexes have correct values
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterator<Node> complex = graphDb.findNodes(Label.label("Complex"));
            while(complex.hasNext()){
                Node cplx = complex.next();
                if((cplx.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("Insulin:p-6Y-Insulin receptor")){
                    String support_score_0 = cplx.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = cplx.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.515, Double.valueOf(support_score_0));
                    assertEquals(4.575, Double.valueOf(abundance_score_0));
                }else if((cplx.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("Insulin receptor")){
                    String support_score_0 = cplx.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = cplx.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.88, Double.valueOf(support_score_0));
                    assertEquals(8.05, Double.valueOf(abundance_score_0));
                }else if((cplx.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("Insulin:Insulin receptor")){
                    String support_score_0 = cplx.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = cplx.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(0.44, Double.valueOf(support_score_0));
                    assertEquals(8.05, Double.valueOf(abundance_score_0));
                }else if((cplx.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("p-6Y-insulin receptor")){
                    String support_score_0 = cplx.getProperty("SUPPORT_SCORE_0").toString();
                    String abundance_score_0 = cplx.getProperty("ABUNDANCE_SCORE_0").toString();
                    assertEquals(1.03, Double.valueOf(support_score_0));
                    assertEquals(4.575, Double.valueOf(abundance_score_0));
                }else if((cplx.getProperty(PropertyType.DISPLAY_NAME.toString()).toString()).equalsIgnoreCase("Insulin")){
                    if(cplx.hasProperty("SUPPORT_SCORE_0")){
                        fail();
                    }else if(cplx.hasProperty("ABUNDANCE_SCORE_0")){
                        fail();
                    }
                }
            }
            tx.success();
        }

        // delete temp directory
        /*
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

         */
    }

    @Test
    void testAddRelWeights_abundance(){

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "HighestSupport"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }


        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterable<Relationship> allRelationships = graphDb.getAllRelationships();
            for(Relationship rel: allRelationships) {
                if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("phospho-p-3-INSR(400-1382)")) { // 3 phosns
                    assertEquals(6.95, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_0").toString()));
                    assertEquals(0.42, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_0").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("phospho-p-2-INSR(28-758)")) { // 2 phosns
                    assertEquals(9.7, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_0").toString()));
                    assertEquals(0.52, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_0").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("INSR(763-1382)")) { // unmodified
                    assertEquals(0.0, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_0").toString()));
                    assertEquals(0.72, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_0").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("4xHC-INS(90-110)")) { // not mapped
                    assertEquals(12.9, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_0").toString()));
                    assertEquals(1.5, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_0").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("Insulin:Insulin receptor")) {
                    assertEquals(8.875, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_0").toString()));
                    assertEquals(1.06, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_0").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("Insulin:p-6Y-Insulin receptor")) {
                    assertEquals(10.612, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_0").toString()));
                    assertEquals(0.985, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_0").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("p-6Y-insulin receptor")) {
                    assertEquals(8.325, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_0").toString()));
                    assertEquals(0.47, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_0").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("Insulin receptor")) {
                    assertEquals(4.85, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_0").toString()));
                    assertEquals(0.62, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_0").toString()));
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
    void testAddRelWeights_lfc(){

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

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {

            InputStream sysInBackup2 = System.in; // backup System.in to restore it later
            ByteArrayInputStream in2 = new ByteArrayInputStream("UniProtID mod_Seq pVal expr".getBytes());
            System.setIn(in2);
            edb.mapMQPhosphopeps(TEST_DATA_FILE, "HighestSupport"); //TODO change to actual input file
            System.setIn(sysInBackup2);// reset System.in to its original

        } catch (IOException e) {
            e.printStackTrace();
        }



        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterable<Relationship> allRelationships = graphDb.getAllRelationships();
            for(Relationship rel: allRelationships) {
                if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("phospho-p-3-INSR(400-1382)")) { // 3 phosns
                    assertEquals(0.69, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_1").toString()));
                    assertEquals(0.42, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_1").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("phospho-p-2-INSR(28-758)")) { // 2 phosns
                    assertEquals(1.12, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_1").toString()));
                    assertEquals(0.52, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_1").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("INSR(763-1382)")) { // unmodified
                    assertEquals(0.0, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_1").toString()));
                    assertEquals(0.72, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_1").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("4xHC-INS(90-110)")) { // not mapped
                    assertEquals(1.2, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_1").toString()));
                    assertEquals(1.5, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_1").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("Insulin:Insulin receptor")) {
                    assertEquals(0.92, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_1").toString()));
                    assertEquals(1.06, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_1").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("Insulin:p-6Y-Insulin receptor")) {
                    assertEquals(1.092, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_1").toString()));
                    assertEquals(0.985, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_1").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("p-6Y-insulin receptor")) {
                    assertEquals(0.985, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_1").toString()));
                    assertEquals(0.47, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_1").toString()));
                } else if (rel.getEndNode().getProperty(PropertyType.DISPLAY_NAME.toString()).equals("Insulin receptor")) {
                    assertEquals(0.64, Double.parseDouble(rel.getProperty("WEIGHT_ABUNDANCE_1").toString()));
                    assertEquals(0.62, Double.parseDouble(rel.getProperty("WEIGHT_SUPPORT_1").toString()));
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
    ///////////////////////////////////////////////////////////////////// DONE /////////////////////////////////////////////////////////////////////

    @Test
    void mapProtein_ScoredBy1() { // when we have a protein with multiple uids and one of them maps but the other doesnt

    }

    @Test
    void mapProtein_ScoredByBoth() { // when we have a protein with multiple uids and both map (take highest)

    }

    @Test
    void pathwayReport(){}


    @Test
    void testWriteSIF(){

        File tempDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraph = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");

        // delete temp graph
        try{
            FileUtils.deleteDirectory(tempDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraph);

        try(Transaction tx = graphDb.beginTx()) {
            Node p1 = graphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            p1.setProperty(PropertyType.DISPLAY_NAME.toString(), "p1");
            p1.setProperty(PropertyType.TYPE.toString(), "Protein");
            p1.setProperty(PropertyType.LOCATION.toString(), "Cytosol");
            p1.setProperty("SUPPORT_SCORE_1", 1);
            p1.setProperty("SUPPORT_SCORE_2", 0.9);
            p1.setProperty("ABUNDANCE_SCORE_1", 0.1);
            p1.setProperty("ABUNDANCE_SCORE_2", 0.2);

            Node p2 = graphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            p2.setProperty(PropertyType.DISPLAY_NAME.toString(), "p2");
            p2.setProperty(PropertyType.TYPE.toString(), "Protein");
            p2.setProperty(PropertyType.LOCATION.toString(), "Cytosol");
            p2.setProperty(PropertyType.INTEGRATED.toString(), "True");
            p2.setProperty("SUPPORT_SCORE_1", 1);
            p2.setProperty("SUPPORT_SCORE_2", 1.5);
            p2.setProperty("ABUNDANCE_SCORE_1", 0.1);
            p2.setProperty("ABUNDANCE_SCORE_2", 0.2);

            Node sm1 = graphDb.createNode(Label.label(LabelTypes.PHYSICAL_ENTITY.toString()));
            sm1.setProperty(PropertyType.DISPLAY_NAME.toString(), "sm1");
            sm1.setProperty(PropertyType.TYPE.toString(), "Small_Molecule");
            sm1.setProperty(PropertyType.LOCATION.toString(), "Cytosol");

            Node phos1 = graphDb.createNode(Label.label(LabelTypes.PHOSPHORYLATION.toString()));
            phos1.setProperty(PropertyType.DISPLAY_NAME.toString(), "phos1");
            phos1.setProperty(PropertyType.TYPE.toString(), "T");
            phos1.setProperty(PropertyType.LOCATION.toString(), "26");

            Node uid1 = graphDb.createNode(Label.label(LabelTypes.UNIPROT_ID.toString()));
            uid1.setProperty(PropertyType.DISPLAY_NAME.toString(), "uid1");

            Node rxn1 = graphDb.createNode(Label.label(LabelTypes.INTERACTION.toString()));
            rxn1.setProperty(PropertyType.DISPLAY_NAME.toString(), "rxn1");
            rxn1.setProperty(PropertyType.TYPE.toString(), "BchmRXN");
            rxn1.setProperty(PropertyType.LOCATION.toString(), "Cytosol");

            Relationship r1 = p1.createRelationshipTo(rxn1, RelTypes.INPUT);
            r1.setProperty("WEIGHT_SUPPORT_1", 1.5);
            r1.setProperty("WEIGHT_ABUNDANCE_1", 0.1);
            r1.setProperty("WEIGHT_SUPPORT_2", 1.5);
            r1.setProperty("WEIGHT_ABUNDANCE_2", 0.2);
            Relationship r2 = sm1.createRelationshipTo(rxn1, RelTypes.INPUT);
            r2.setProperty("WEIGHT_SUPPORT_1", 1.5);
            r2.setProperty("WEIGHT_ABUNDANCE_1", 0.1);
            r2.setProperty("WEIGHT_SUPPORT_2", 1.5);
            r2.setProperty("WEIGHT_ABUNDANCE_2", 0.2);
            Relationship r3 = rxn1.createRelationshipTo(p2, RelTypes.OUTPUT);
            r3.setProperty("WEIGHT_SUPPORT_1", 1);
            r3.setProperty("WEIGHT_ABUNDANCE_1", 0.1);
            r3.setProperty("WEIGHT_SUPPORT_2", 1.5);
            r3.setProperty("WEIGHT_ABUNDANCE_2", 0.2);
            Relationship r4 = phos1.createRelationshipTo(p2, RelTypes.PHOSPHORYLATION);
            r4.setProperty("WEIGHT_SUPPORT_1", 1.5);
            r4.setProperty("WEIGHT_ABUNDANCE_1", 0.1);
            r4.setProperty("WEIGHT_SUPPORT_2", 1.5);
            r4.setProperty("WEIGHT_ABUNDANCE_2", 0.2);
            Relationship r5 = uid1.createRelationshipTo(p1, RelTypes.ID_BELONGS_TO);
            r5.setProperty("WEIGHT_SUPPORT_1", 1.5);
            r5.setProperty("WEIGHT_ABUNDANCE_1", 0.5);
            r5.setProperty("WEIGHT_SUPPORT_2", 1.5);
            r5.setProperty("WEIGHT_ABUNDANCE_2", 0.5);
            Relationship r6 = uid1.createRelationshipTo(p2, RelTypes.ID_BELONGS_TO);
            r6.setProperty("WEIGHT_SUPPORT_1", 1.5);
            r6.setProperty("WEIGHT_ABUNDANCE_1", 0.5);
            r6.setProperty("WEIGHT_SUPPORT_2", 1.5);
            r6.setProperty("WEIGHT_ABUNDANCE_2", 0.5);

            for (String allPropertyKey : graphDb.getAllPropertyKeys()) {
                System.out.println(allPropertyKey);
            }
            tx.success();
        }
        graphDb.shutdown();

        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraph, tempDir);

        try {
            edb.writeSIF();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // TODO compare printed file to file from expected file dir


        // delete temp graph
//        try{
//            FileUtils.deleteDirectory(tempDir);
//        }catch (IOException ex){
//            ex.printStackTrace();
//        }
    } // print to file then use the file comparison thing

    @Test
    void test(){
        File tempOutputDir = new File("/Users/huckstep.h/Documents/neo4j/PhlashyName_home/temp/");
        File tempGraphDir = new File("/Users/huckstep.h/Documents/neo4j/PhlashyName_home/temp/HUMAN/");
        File seansData = new File("/Users/huckstep.h/Documents/neo4j/PhlashyName_home/tutorial/SEANSDATA.tsv");

        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            ResourceIterable<Label> allLabels = graphDb.getAllLabels();
            for (Label lbl: allLabels) {
                System.out.println(lbl);
            }

            HashSet<String> hs = new HashSet<String>();
            ResourceIterator<Node> smallMolecule = graphDb.findNodes(Label.label("SmallMolecule"));
            for (ResourceIterator<Node> it = smallMolecule; it.hasNext(); ) {
                Node sm = it.next();
                hs.add(sm.getProperty(PropertyType.DISPLAY_NAME.toString()).toString());
            }

            System.out.println(hs.size());


            ResourceIterator<Node> uid = graphDb.findNodes(Label.label("UNIPROT_ID"));
            Integer count = 0;
            for (ResourceIterator<Node> it = uid; it.hasNext(); ) {
                Node id = it.next();
                count ++;
            }

            System.out.println(count);

            System.out.println(hs);
        }
    }

    @Test
    void test2(){ // just printing the db
        File tempOutputDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraphDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");

        // in case graph wasn't already deleted
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        // Make a new graph
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraphDir, true, "human");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original

        // Map data onto it
        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);
        try {
            edb.writeSIF();
        } catch (IOException e) {
            e.printStackTrace();
        }




    }

    @Test
    void testKinaseReport() {

        File tempOutputDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraphDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");

        /*
        // in case graph wasn't already deleted
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        // Make a new graph
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraphDir, false, "human");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original

         */

        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir, tempOutputDir);

        edb.printDatabase();
        //if one line in the report == kinases 1 then we gucci

        // delete temp directory
        /*
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
         */
    }

}