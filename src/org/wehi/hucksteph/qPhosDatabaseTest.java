package org.wehi.hucksteph;

import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.Test;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import static org.junit.jupiter.api.Assertions.*;

class qPhosDatabaseTest {

    private final File TEST_OWL_FILE = new File("/Users/huckstep.h/Dropbox/IdeaProjects/ReactoSitePlus/test/InsulinTestNtwk.owl");

    private final File TEST_DATA_FILE = new File("/Users/huckstep.h/Dropbox/IdeaProjects/ReactoSitePlus/test/mapTestData.txt");
    private final File DATABASE_EXPECTED_PATH = new File("/Users/huckstep.h/Dropbox/IdeaProjects/ReactoSitePlus/test/expected/");
    private final File DATABASE_ACTUAL_PATH = new File("/Users/huckstep.h/Dropbox/IdeaProjects/ReactoSitePlus/test/actual/");

    @Test
    void allDS() {

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

        qPhosDatabase qPhos = new qPhosDatabase(tempGraphDir,tempOutputDir);
        try {
            qPhos.allDS();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    @Test
    void allMappings() {
    }

    @Test
    void allNbhds() {
    }
}