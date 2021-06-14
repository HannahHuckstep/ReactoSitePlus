package org.wehi.hucksteph;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.junit.jupiter.api.Test;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.unsafe.impl.batchimport.input.InputException;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class MeasuredDatabaseTest {

    private final File TEST_OWL_FILE = new File("test/InsulinTestNtwk.owl");
    private final File TEST_DATA_FILE = new File("test/mapTestData.txt");
    private final File DATABASE_EXPECTED_PATH = new File("test/expected/");
    private final File DATABASE_ACTUAL_PATH = new File("test/actual/");

    @Test
    void distributionArrayTest() {

        List<Integer> x = new ArrayList<>(Arrays.asList(100, 1000, 2000, 7000, 10000, 82149, 87324, 32949, 92047, 232456));

        for (Integer integer :x) {

            Integer subsetSize = integer;
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

            assertTrue(subsetSize <= distribution.size());
        }

    }

    @Test
    void traversalTest_DS_UID() {
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.traversal("P06213", "downstream", "0");
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file1 = new File(tempOutputDir + "/P06213_downstream.tsv");
        File file2 = new File("test/expected/traversal/P06213_downstream.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file1, file2));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file3 = new File(tempOutputDir + "/TraversalReport_downstream_P06213.tsv");
        File file4 = new File("test/expected/traversal/TraversalReport_downstream_P06213.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file3, file4));
        } catch (IOException e) {
            e.printStackTrace();
        }

        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
    }

    @Test
    void traversalTest_US_UID() {
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.traversal("Q99943", "upstream", "0");
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file1 = new File(tempOutputDir + "/Q99943_upstream.tsv");
        File file2 = new File("test/expected/traversal/Q99943_upstream.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file1, file2));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file3 = new File(tempOutputDir + "/TraversalReport_upstream_Q99943.tsv");
        File file4 = new File("test/expected/traversal/TraversalReport_upstream_Q99943.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file3, file4));
        } catch (IOException e) {
            e.printStackTrace();
        }


        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }


    }

    @Test
    void traversalTest_DS_nodeID() {
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.traversal("17", "downstream", "0");
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file1 = new File(tempOutputDir + "/17_Downstream.tsv");
        File file2 = new File("test/expected/traversal/17_Downstream.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file1, file2));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file3 = new File(tempOutputDir + "/TraversalReport_downstream_17.tsv");
        File file4 = new File("test/expected/traversal/TraversalReport_downstream_17.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file3, file4));
        } catch (IOException e) {
            e.printStackTrace();
        }


        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    @Test
    void traversalTest_US_nodeID() {
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.traversal("1", "upstream", "0");
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file1 = new File(tempOutputDir + "/1_Upstream.tsv");
        File file2 = new File("test/expected/traversal/1_Upstream.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file1, file2));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file3 = new File(tempOutputDir + "/TraversalReport_upstream_1.tsv");
        File file4 = new File("test/expected/traversal/TraversalReport_upstream_1.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file3, file4));
        } catch (IOException e) {
            e.printStackTrace();
        }


        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
    }

    @Test
    void traversalThrowErrorWhenWrongInput() {
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);

        Exception exception = assertThrows(IllegalArgumentException.class, () ->
        { mdb.traversal("P06213", "stream", "0");
        });
        System.setIn(sysInBackup);// reset System.in to its original


        String expectedMessage = "direction must equal 'upstream' or 'downstream'";
        String actualMessage = exception.getMessage();

        assertTrue(actualMessage.contains(expectedMessage));




        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    @Test
    void getNeighbours() {
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


        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir);
        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(tempGraphDir);
        try(Transaction tx = graphDb.beginTx()){
            HashSet<Node> nodes = mdb.neighbourTraversal(graphDb.getNodeById(16), 4, graphDb);
            HashSet<Node> nodesExpected = new HashSet<>();
            nodesExpected.add(graphDb.getNodeById(16));
            assertEquals(nodesExpected, nodes);

            HashSet<Node> nodes26 = mdb.neighbourTraversal(graphDb.getNodeById(26), 4, graphDb);
            HashSet<Node> nodesExpected26 = new HashSet<>();
            for (int i = 20; i < 29; i++) {
                nodesExpected26.add(graphDb.getNodeById(i));
            }
            assertEquals(nodesExpected26, nodes26);

            HashSet<Node> nodes34 = mdb.neighbourTraversal(graphDb.getNodeById(34), 4, graphDb);
            HashSet<Node> nodesExpected34 = new HashSet<>();
            nodesExpected34.add(graphDb.getNodeById(1));
            nodesExpected34.add(graphDb.getNodeById(29));
            nodesExpected34.add(graphDb.getNodeById(30));
            nodesExpected34.add(graphDb.getNodeById(32));
            nodesExpected34.add(graphDb.getNodeById(34));
            nodesExpected34.add(graphDb.getNodeById(44));
            nodesExpected34.add(graphDb.getNodeById(45));
            assertEquals(nodesExpected34, nodes34);

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
    void SP_UID_UID(){
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.shortestPath("Q9UH92-2", "Q99943", "Support");
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file1 = new File(tempOutputDir + "/Q9UH92-2_to_Q99943_downstream.tsv");
        File file2 = new File("test/expected/shortestPath/Q9UH92-2_to_Q99943_downstream.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file1, file2));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file3 = new File(tempOutputDir + "/ShortestPath_Q9UH92-2_to_Q99943.tsv");
        File file4 = new File("test/expected/shortestPath/ShortestPath_Q9UH92-2_to_Q99943.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file3, file4));
        } catch (IOException e) {
            e.printStackTrace();
        }

        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    @Test
    void SP_nodeID_nodeID(){
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.shortestPath("13", "26", "Abundance");
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file1 = new File(tempOutputDir + "/13_to_26_downstream.tsv");
        File file2 = new File("test/expected/shortestPath/13_to_26_downstream.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file1, file2));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file3 = new File(tempOutputDir + "/ShortestPath_13_to_26.tsv");
        File file4 = new File("test/expected/shortestPath/ShortestPath_13_to_26.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file3, file4));
        } catch (IOException e) {
            e.printStackTrace();
        }

        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    @Test
    void SP_UID_nodeID(){
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.shortestPath("P06213", "26", "Support");
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file1 = new File(tempOutputDir + "/P06213_to_26_downstream.tsv");
        File file2 = new File("test/expected/shortestPath/P06213_to_26_downstream.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file1, file2));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file3 = new File(tempOutputDir + "/ShortestPath_P06213_to_26.tsv");
        File file4 = new File("test/expected/shortestPath/ShortestPath_P06213_to_26.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file3, file4));
        } catch (IOException e) {
            e.printStackTrace();
        }

        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    @Test
    void SP_SupportScore(){
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.shortestPath("P06213", "22", "Support");
        } catch (IOException e) {
            e.printStackTrace();
        }
        File file1 = new File(tempOutputDir + "/P06213_to_22_downstream.tsv");
        File file2 = new File("test/expected/shortestPath/P06213_to_22_downstream_SS.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file1, file2));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file3 = new File(tempOutputDir + "/ShortestPath_P06213_to_22.tsv");
        File file4 = new File("test/expected/shortestPath/ShortestPath_P06213_to_22_SS.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file3, file4));
        } catch (IOException e) {
            e.printStackTrace();
        }

        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    @Test
    void SP_AbundanceScore(){
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.shortestPath("P06213", "22", "Abundance");
        } catch (IOException e) {
            e.printStackTrace();
        }
        File file1 = new File(tempOutputDir + "/P06213_to_22_downstream.tsv");
        File file2 = new File("test/expected/shortestPath/P06213_to_22_downstream.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file1, file2));
        } catch (IOException e) {
            e.printStackTrace();
        }

        File file3 = new File(tempOutputDir + "/ShortestPath_P06213_to_22.tsv");
        File file4 = new File("test/expected/shortestPath/ShortestPath_P06213_to_22.tsv");

        try {
            assertTrue(FileUtils.contentEquals(file3, file4));
        } catch (IOException e) {
            e.printStackTrace();
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
    void SP_nodeID_UID_noPath(){
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);

        Exception exception = assertThrows(Exception.class, () ->
        { mdb.shortestPath("34", "Q99943", "Support");
        });
        System.setIn(sysInBackup);// reset System.in to its original

        String expectedMessage = "No path between nodes 34 and Q99943 either upstream or downstream";
        String actualMessage = exception.getMessage();

        assertTrue(actualMessage.contains(expectedMessage));


        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    @Test
    void SP_throwErrorWhenUIDNotInDb(){
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);

        Exception exception = assertThrows(Exception.class, () ->
        { mdb.shortestPath("P12345", "Q99943", "Support");
        });
        System.setIn(sysInBackup);// reset System.in to its original


        String expectedMessage = "No path between nodes 34 and Q99943 either upstream or downstream";
        String actualMessage = exception.getMessage();

        assertTrue(actualMessage.contains(expectedMessage));


        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
    }

    @Test
    void SP_throwErrorWhenNodeIDNotInDb(){
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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);

        Exception exception = assertThrows(Exception.class, () ->
        { mdb.shortestPath("100", "Q99943", "Support");
        });
        System.setIn(sysInBackup);// reset System.in to its original


        String expectedMessage = "No path between nodes 34 and Q99943 either upstream or downstream";
        String actualMessage = exception.getMessage();

        assertTrue(actualMessage.contains(expectedMessage));


        // delete temp directory
        try{
            FileUtils.deleteDirectory(tempOutputDir);
        }catch (IOException ex){
            ex.printStackTrace();
        }
    }

    @Test
    void resetScores(){

    }

    @Test
    void TraverseNonProtNodeID(){

    }

    @Test
    void nbhdStats(){
        File tempOutputDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraphDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");
        File qPhosFile = new File("");

        // in case graph wasn't already deleted
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        // Make a new graph
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraphDir, false, "human");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original


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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);

        try {
            mdb.empiricalNullDistribution(qPhosFile, 4,10, 10);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    @Test
    void nbhd(){
        File tempOutputDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/");
        File tempGraphDir = new File(DATABASE_ACTUAL_PATH+ "/toBeDeleted/GRAPH/");
        File qPhosFile = new File("");

        // in case graph wasn't already deleted
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        // Make a new graph
        DatabaseFactory dbf = new DatabaseFactory(TEST_OWL_FILE, tempGraphDir, false, "human");
        dbf.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original


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

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);

        try {
            mdb.nbhdAnalysis(4, "0");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    @Test
    void MCNTest(){

    }

    @Test
    void myDist(){
        File tempOutputDir = new File("/Users/huckstep.h/Documents/neo4j/PhlashyName/liveDemo/");
        File tempGraphDir = new File( "/Users/huckstep.h/Documents/neo4j/PhlashyName/liveDemo/HUMAN/");
        File miniDataFile = new File( "/Users/huckstep.h/Documents/neo4j/Project2/empiricalDist/qPhos_all_data_smaller.txt");

        MeasuredDatabase mdb = new MeasuredDatabase(tempGraphDir, tempOutputDir);
        try {
            mdb.empiricalNullDistribution(miniDataFile, 4, 1000, 3);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    @Test
    void test(){
        File tempOutputDir = new File("test/actual/");
        //File tempGraphDir = new File( "test/actual/toBeDeleted/GRAPH/");
        File tempGraphDir = new File( "/Users/huckstep.h/Documents/neo4j/PhlashyName/liveDemo/HUMAN/");

        EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(tempGraphDir,tempOutputDir );

        try {
            edb.writeSIF();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}