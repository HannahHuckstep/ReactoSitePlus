package org.wehi.hucksteph;

import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.Test;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import static org.junit.jupiter.api.Assertions.*;

class ReactomeDatabaseTest {

    File INTEGRATION_RXM_OWL = new File("test/integrationRXM.owl");
    File INTEGRATION_PSP_OWL = new File("test/integrationPSP.owl");
    File EXPECTED_INTEGRATION = new File("test/expected/integrated/");
    File ACTUAL_OUTPUT = new File("test/actual/integrated/");

    @Test
    /*
    Tests performed in this test case:
    - if there is a phosphorylated protein in RXM not involved in a bchm rxn, should make one
    - The UID controlling a rxn in PSP is not in RXM, should create a new UID, protein and controller node
    - A phosphorylation in PSP is not in RXM, should create a new prot node, with new phos, the reaction node & the appropriate controls
    - A phosphorylated protein in PSP is not in RXM, should create new UID, prot, phos, and rxn. Then attach to appropriate controls
    - A rxn in RXM doesnt have a catalysis controlling it, should create a catalysis reaction and add the appropriate controller
    - The RXM control does not match the PSP control, should create a new separate reaction, phosphorylations, catalysis and control
    - No appropriate 'before' (unphosphrylated) protein in the same cellular location as the rxn, should make a new one and link to appropriate UID
    - A controller exists for a reaction but not in the right cellular location, should create one ans link to appropriate UID

    All of these cases have been built into the Reactome (RXM) and PhosphoSitePlus (PSP) owl files
    Each is addressed in the answer integrated Neo4j Network (ACTUAL_OUTPUT)
     */
    public void basicIntegratePSP3() throws IOException {
        File ACTUAL_RXM_GRAPH = new File(ACTUAL_OUTPUT + "/RXMGRAPH/");
        File ACTUAL_PSP_GRAPH = new File(ACTUAL_OUTPUT + "/PSPGRAPH/");

        //Reactome Graph
        InputStream sysInBackup = System.in; // backup System.in to restore it later
        ByteArrayInputStream in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        DatabaseFactory rdb = new DatabaseFactory(INTEGRATION_RXM_OWL,ACTUAL_RXM_GRAPH, false, "h");
        rdb.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original

        //PhosphoSitePlus Graph
        sysInBackup = System.in; // backup System.in to restore it later
        in = new ByteArrayInputStream("y".getBytes());
        System.setIn(in);
        DatabaseFactory pdb = new DatabaseFactory(INTEGRATION_PSP_OWL,ACTUAL_PSP_GRAPH, false, "h");
        pdb.createDBfromOWL();
        System.setIn(sysInBackup);// reset System.in to its original

        // just made
        ReactomeDatabase rdb_integrated = new ReactomeDatabase(ACTUAL_RXM_GRAPH, ACTUAL_OUTPUT);
        rdb_integrated.IntegratePSP(ACTUAL_PSP_GRAPH);

        // just made
        GraphDatabaseService actual = new GraphDatabaseFactory().newEmbeddedDatabase(ACTUAL_RXM_GRAPH);
        // answer
        File expectedGraph = new File(EXPECTED_INTEGRATION + "/integratedGRAPH/");
        GraphDatabaseService expected = new GraphDatabaseFactory().newEmbeddedDatabase(expectedGraph);

        // cypher query that returns everything
        String a = null;
        try(Transaction tx = actual.beginTx()){
            Result result = actual.execute("MATCH (n) RETURN n");
            a = result.resultAsString();
            tx.success();
        }
        String e = null;
        try(Transaction tx = expected.beginTx()){
            Result result = expected.execute("MATCH (n) RETURN n");
            e = result.resultAsString();
            tx.success();
        }

        assertEquals(e,a);


        File actualIntegrationStats = new File(ACTUAL_OUTPUT + "/IntegrationStats.txt");
        File expectedIntegrationStats = new File(EXPECTED_INTEGRATION + "/IntegrationStats.txt");
        boolean b = FileUtils.contentEquals(expectedIntegrationStats, actualIntegrationStats);
        if(!b){
            fail();
        }

        try{
            FileUtils.deleteDirectory(ACTUAL_OUTPUT);
        }catch (IOException ex){
            ex.printStackTrace();
        }

    }

    
}