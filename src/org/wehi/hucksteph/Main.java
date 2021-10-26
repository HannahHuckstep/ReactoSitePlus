package org.wehi.hucksteph;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import org.neo4j.unsafe.impl.batchimport.input.InputException;
import scala.Int;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Main {

    public static void main(String[] args) throws IOException {




        ArgumentParser parser = ArgumentParsers.newFor("org.wehi.hucksteph").build().defaultHelp(true)
                .description("A tool for creating, modifying and interrogating phosphoproteomics networks");

        parser.addArgument("--mode", "-m")
                .dest("mode")
                .help("\nThe function you would like to perform. \n" +
                        "Options are:\n" +
                                "\"PrintDatabase\", takes in input database [-idb]\n" +
                                "\"AmountWithLabel\", takes in input database [-idb] and the name of the label of interest [-l]\n" +
                                "\"WriteAllUIDs\", takes an input database [-idb] and an output path [-op]\n" +
                                "\"WritePhos\", takes an input database [-idb] and an output path [-op]\n" +
                                "\"WriteDBtoSIF\", takes an input database [-idb] and an output path [-op]\n" +
                                "\"IntegratePSP\", takes in input Reactome database [-idb], takes the PSP database [-psp], an output path [-op], and a species [-s]\n" +
                                "\"BinomialNeighbourhoodAnalysis\", takes in a measured input database [-idb], an output path [-op], and the depth of the traversal [-d]\n"+
                                "\"MapPeptides\", takes an input database [-idb] and an output path [-op], a file to map onto the database [-idf], and the optional Abundance Score mapping method preferred [-as] (\"HighestSupport\" is defalut)\n" +
                                "\"CreateDB\", takes an OWL file [-iof], an output path [-op], an optional update boolean [-u] (can be T or F, default is T), and the species of graph you'd like to make [-s] (can be human (h) or mouse(m))\n" +
                                "\"NeighbourhoodAnalysis\", takes in a measured input database [-idb], an output path [-op], the depth of the traversal [-d], the experiment name of interest [-en], and the file containing the pre-calculated empirical distribution per neighbourhood [-idf]\n"+
                                "\"ShortestPath\", takes in a measured input database [-idb], an output path [-op], a starting node id [-sid], a ending node id [-eid], and the weight type to be traversed [-ew] (can be either \"Abundance\" (a) or \"Support\" (s))\n"+
                                "\"MinimalConnectionNetwork\", takes in a measured input database [-idb], an output path [-op], and the experiment name of interest [-en]\n"+
                                "\"TraversalAnalysis\", takes in a measured input database [-idb], an output path [-op], a UniProt ID or database ID to look downstream of [-p], the direction of the traversal [-dir], and the experiment name of interest [-en]\n" +
                                "\nqPhosDs"+
                                "\nqPhosMapALL\t [idb][op][idf]" +
                                "\nqPhosMap \t [idb][op][idf][en]" +
                                "\nqPhosNbhd" +
                                "\nqPhosMCN \t [idb][op][en]" +
                                "\nqPhosED \t [idb][op][idf][d][ss][rn]"+
                                "\nmouseED \t [idb][op][idf][d][ss][rn]"+
                        "\n"

                        )
                .type(String.class)
                .choices("CreateDB", //done
                        "PrintDatabase", //done
                        "WriteAllUIDs", //done
                        "WritePhos", //done
                        "WriteDBtoSIF", //done
                        "MapPeptides", //done
                        "IntegratePSP", //done
                        "BinomialNeighbourhoodAnalysis", //done
                        "NeighbourhoodAnalysis",
                        "TraversalAnalysis", //done
                        "ShortestPath", //done
                        "MinimalConnectionNetwork", //done
                        "ResetScores",
                        "GetSpecies",
                        "AmountWithLabel", //done
                        "PrintAllProperties",
                        "qPhosDs",
                        "qPhosMap",
                        "qPhosNbhd",
                        "qPhosMCN",
                        "qPhosED",
                        "mouseED"

                )
                .required(true);
        parser.addArgument("--input_owl_file", "-iof")
                .dest("input_owl_file")
                .nargs("?")
                .help("The OWL file to input")
                .type(File.class);
        parser.addArgument("--output_path", "-op")
                .dest("output_path")
                .nargs("?")
                .help("The destination path of the file");
        parser.addArgument("--update", "-u")
                .dest("update")
                .nargs("?")
                .help("For CreatDB use flag if you would like the database to be updated while being created");
        parser.addArgument("--species", "-s" )
                .dest("species")
                .nargs("?")
                .help("The Species of netowrk to build -human (h) or mouse (m)");
        parser.addArgument("--input_db", "-idb")
                .dest("input_db")
                .nargs("?")
                .help("The graph database directory");
        parser.addArgument("--input_data_file", "-idf")
                .dest("input_data_file")
                .nargs("?")
                .help("The input data file to map");
        parser.addArgument("--input_data_file2", "-idf2")
                .dest("input_data_file2")
                .nargs("?")
                .help("The second input data file to map");
        parser.addArgument("--experiment_name", "-en" )
                .dest("experiment")
                .nargs("?")
                .help("The experiment name to use for the  analysis");
        parser.addArgument("--abundanceScore", "-as" )
                .dest("abundanceScore")
                .nargs("?")
                .help("The abundance score mapping method preferred");
        parser.addArgument("--psp_db", "-psp" )
                .dest("psp")
                .nargs("?")
                .help("The directory holding the PhosphositePlus neo4j Database to integrate");
        parser.addArgument("--depth", "-d" )
                .dest("depth")
                .nargs("?")
                .help("The depth of the neighbourhood traversal");
        parser.addArgument("--direction", "-dir" )
                .dest("direction")
                .nargs("?")
                .help("The direction of the traversal");
        parser.addArgument("--protein", "-p" )
                .dest("protein")
                .nargs("?")
                .help("The protein to start the traversal");
        parser.addArgument("--startID", "-sid" )
                .dest("startID")
                .nargs("?")
                .help("The id of the starting node to get the shortest path");
        parser.addArgument("--endID", "-eid" )
                .dest("endID")
                .nargs("?")
                .help("The id of the ending node to get the shortest path");
        parser.addArgument("--label", "-l" )
                .dest("label")
                .nargs("?")
                .help("The name of the label of interest");
        parser.addArgument("--edge_weights", "-ew" )
                .dest("edge_weights")
                .nargs("?")
                .help("The edge weights to be traversed");
        parser.addArgument("--subsetSize", "-ss" )
                .dest("subset_size")
                .nargs("?")
                .help("The # of proteins to be mapped");
        parser.addArgument("--repetitionNum", "-rn" )
                .dest("repetition_num")
                .nargs("?")
                .help("The number of times the experiment should be repeated");

        try{

            Namespace ns = parser.parseArgs(args);
            String mode = ns.get("mode");

            if(mode.equalsIgnoreCase("CreateDB")){
                if(ns.getAttrs().get("input_owl_file") == null){
                    throw new NullPointerException("Missing the input OWL file");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                } else if(ns.getAttrs().get("species") == null){
                    throw new NullPointerException("Missing species specified");
                    // if update param is missing default is true
                }else if(ns.getAttrs().get("update") == null){
                    File input_owl_file = new File(ns.get("input_owl_file").toString());
                    File output_graph_directory = new File(ns.get("output_path").toString());
                    DatabaseFactory db = new DatabaseFactory(input_owl_file,output_graph_directory);
                    db.createDBfromOWL();

                }else{
                    File input_owl_file = new File(ns.get("input_owl_file").toString());
                    File output_graph_directory = new File(ns.get("output_path").toString());
                    String update = ns.get("update");
                    String species = ns.get("species");
                    if (update.equalsIgnoreCase("True") | update.equalsIgnoreCase("T")){
                        DatabaseFactory db = new DatabaseFactory(input_owl_file, output_graph_directory, true, species);
                        db.createDBfromOWL();
                    }else if (update.equalsIgnoreCase("false") | update.equalsIgnoreCase("F")){
                        DatabaseFactory db = new DatabaseFactory(input_owl_file, output_graph_directory, false, species);
                        db.createDBfromOWL();
                    }else{
                        throw new IllegalArgumentException("Update parameter must be True or False");
                    }
                }
            }
            else if(mode.equalsIgnoreCase("AmountWithLabel")) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                }else if(ns.getAttrs().get("label") == null){
                    throw new NullPointerException("Missing Label of interest");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    String label = ns.get("label");
                    EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(input_db);
                    edb.printLabelNumber(label);
                }
            }
            else if(mode.equalsIgnoreCase("PrintDatabase")) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(input_db);
                    edb.printDatabase();
                }
            }
            else if(mode.equalsIgnoreCase("WriteAllUIDs")) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(input_db, output_path);
                    try {
                        edb.printAllUniProtIDs();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
            else if(mode.equalsIgnoreCase("WritePhos")) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(input_db, output_path);
                    try {
                        edb.printAllPhosphorylations();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
            else if(mode.equalsIgnoreCase("WriteDBtoSIF")) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    EmbeddedNeo4jDatabase edb = new EmbeddedNeo4jDatabase(input_db, output_path);
                    try {
                        edb.writeSIF();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
            else if(mode.equalsIgnoreCase("IntegratePSP")) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                } else if(ns.getAttrs().get("psp") == null){
                    throw new NullPointerException("Missing PhosphositePlus Database to integrate");
                } else if(ns.getAttrs().get("species") == null){
                    throw new NullPointerException("Missing species specified");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    File psp = new File(ns.get("psp").toString());
                    String species = ns.get("species");
                    ReactomeDatabase rxmdb = new ReactomeDatabase(input_db, output_path);
                    rxmdb.IntegratePSP(psp, species);
                }
            }
            else if(mode.equalsIgnoreCase("MapPeptides")) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                } else if(ns.getAttrs().get("input_data_file") == null){
                    throw new NullPointerException("Missing Data file to map");
                }else if(ns.getAttrs().get("abundanceScore") == null){
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    File p2p = new File(ns.get("input_data_file").toString());
                    MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                    mdb.mapMQPhosphopeps(p2p, "HighestSupport");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    File p2p = new File(ns.get("input_data_file").toString());
                    String method = ns.getString("abundanceScore");
                    MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                    mdb.mapMQPhosphopeps(p2p, method);
                }
            }
            else if(mode.equalsIgnoreCase("BinomialNeighbourhoodAnalysis")) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                } else if(ns.getAttrs().get("depth") == null) {
                    throw new NullPointerException("Missing depth parameter for neighbourhood traversal");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    String depth = ns.getString("depth");
                    Integer d = 0;
                    try{
                         d = Integer.valueOf(depth);
                    }catch (NumberFormatException e){
                        System.out.println("Depth must be an integer");
                    }

                    MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                    mdb.binomialNeighbourhood(d);
                }
            }
            else if(mode.equalsIgnoreCase("NeighbourhoodAnalysis")) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                } else if(ns.getAttrs().get("depth") == null) {
                    throw new NullPointerException("Missing depth parameter for neighbourhood traversal");
                }else if(ns.getAttrs().get("experiment") == null) {
                    throw new NullPointerException("Missing experiment parameter for network traversal");
                }else if(ns.getAttrs().get("input_data_file") == null){
                    throw new NullPointerException("Missing Data file to map");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    String expt = ns.getString("experiment");
                    String depth = ns.getString("depth");
                    Integer d = 0;
                    File ed = new File(ns.get("input_data_file").toString());
                    try{
                        d = Integer.valueOf(depth);
                    }catch (NumberFormatException e){
                        System.out.println("Depth must be an integer");
                        System.exit(1);
                    }
                    if(ns.getAttrs().get("input_data_file2") == null){
                        MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                        mdb.nbhdAnalysis(d, expt, ed);
                    }else{
                        File ed2 = new File(ns.get("input_data_file2").toString());
                        MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                        mdb.nbhdAnalysis(d, expt, ed, ed2);
                    }

                }
            }
            else if(mode.equalsIgnoreCase("TraversalAnalysis") ) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                } else if(ns.getAttrs().get("direction") == null) {
                    throw new NullPointerException("Missing direction parameter for network traversal");
                } else if(ns.getAttrs().get("protein") == null) {
                    throw new NullPointerException("Missing protein parameter for network traversal");
                } else if(ns.getAttrs().get("experiment") == null) {
                    throw new NullPointerException("Missing experiment parameter for network traversal");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    String direction = ns.getString("direction");
                    String protein = ns.getString("protein");
                    String experiment = ns.getString("experiment");
                    if(!direction.equalsIgnoreCase("upstream") & !direction.equalsIgnoreCase("downstream")) {
                        throw new InputException("Direction must be either \"upstream\" or \"downstream\" ");
                    }
                    MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                    mdb.traversal(protein, direction, experiment);
                }
            }
            else if(mode.equalsIgnoreCase("ShortestPath")  ) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                } else if(ns.getAttrs().get("startID") == null) {
                    throw new NullPointerException("Missing node ID for the start of the shortest path");
                } else if(ns.getAttrs().get("endID") == null) {
                throw new NullPointerException("Missing node ID for the end of the shortest path");
                } else if(ns.getAttrs().get("edge_weights") == null){
                    throw new NullPointerException("Missing the type of weights to be traversed");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    String s = ns.getString("startID");
                    String e = ns.getString("endID");
                    String weight = ns.get("edge_weights");
                    if(weight.equalsIgnoreCase("Abundance")| weight.equalsIgnoreCase("Support") | weight.equalsIgnoreCase("a")| weight.equalsIgnoreCase("s")){
                        MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                        mdb.shortestPath(s,e, weight);
                    }else{
                        throw new InputException("Weights to be traversed can be either \"Abundance\" (\"a\")or \"Support\" (\"s\")");
                    }

                }
            }
            else if(mode.equalsIgnoreCase("MinimalConnectionNetwork")  ) {
                if (ns.getAttrs().get("input_db") == null) {
                    throw new NullPointerException("Missing the input: database directory");
                } else if(ns.getAttrs().get("output_path") == null){
                    throw new NullPointerException("Missing output path to write to");
                //} else if(ns.getAttrs().get("edge_weights") == null){
                //    throw new NullPointerException("Missing the type of weights to be traversed");
                } else if(ns.getAttrs().get("experiment") == null) {
                    throw new NullPointerException("Missing experiment parameter for network traversal");
                }else{
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                    String weight = ns.get("edge_weights");
                    String e = ns.get("experiment");
                    //if(weight.equalsIgnoreCase("Abundance")| weight.equalsIgnoreCase("Support") | weight.equalsIgnoreCase("a")| weight.equalsIgnoreCase("s")){
                        MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                        mdb.minimalConnectionNetwork(e);
                    //}ese{
                    //    throw new InputException("Weights to be traversed can be either \"Abundance\" (\"a\")or \"Support\" (\"s\")");
                    //}

                }
            }
            /////////////////////////////////////////////////////////////////////////////////
            else if(mode.equalsIgnoreCase("qPhosDs")  ) {
                    File input_db = new File(ns.get("input_db").toString());
                    File output_path = new File(ns.get("output_path").toString());
                qPhosDatabase qPhos = new qPhosDatabase(input_db,output_path);
                qPhos.allDS();
            }
            else if(mode.equalsIgnoreCase("qPhosMapAll")  ) {
                File input_db = new File(ns.get("input_db").toString());
                File output_path = new File(ns.get("output_path").toString());
                File p2p = new File(ns.get("input_data_file").toString());
                qPhosDatabase qPhos = new qPhosDatabase(input_db,output_path);
                qPhos.allMappings(p2p);
            }
            else if(mode.equalsIgnoreCase("qPhosMap")  ) {
                File input_db = new File(ns.get("input_db").toString());
                File output_path = new File(ns.get("output_path").toString());
                File p2p = new File(ns.get("input_data_file").toString());
                String expr = ns.getString("experiment");
                qPhosDatabase qPhos = new qPhosDatabase(input_db,output_path);
                qPhos.allMappings(p2p, expr);
            }
            else if(mode.equalsIgnoreCase("qPhosNbhd")  ) {
                File input_db = new File(ns.get("input_db").toString());
                File output_path = new File(ns.get("output_path").toString());
                qPhosDatabase qPhos = new qPhosDatabase(input_db,output_path);
                qPhos.allNbhds(4);

            }
            else if(mode.equalsIgnoreCase("qPhosMCN")  ) {
                File input_db = new File(ns.get("input_db").toString());
                File output_path = new File(ns.get("output_path").toString());
                String experiment = ns.getString("experiment");
                qPhosDatabase qPhos = new qPhosDatabase(input_db,output_path);
                qPhos.qPhosMinimalConnectionNetwork(experiment);

            }
            else if(mode.equalsIgnoreCase("qPhosED")  ) { //empirical distribution generator
                // [idb][op][idf][d][ss][rn]
                File input_db = new File(ns.get("input_db").toString());
                File output_path = new File(ns.get("output_path").toString());
                File data_file = new File(ns.get("input_data_file").toString());
                String subsetSize = ns.getString("subset_size");
                String repetitionNum = ns.getString("repetition_num");
                String depth = ns.getString("depth");
                Integer d = 0;
                try{
                    d = Integer.valueOf(depth);
                }catch (NumberFormatException e){
                    System.out.println("Depth must be an integer");
                }
                Integer ss = Integer.valueOf(subsetSize);
                Integer rn = Integer.valueOf(repetitionNum);
                MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                mdb.empiricalNullDistribution(data_file,d, ss, rn);
            }
            else if(mode.equalsIgnoreCase("mouseED")  ) { //empirical distribution generator
                // [idb][op][idf][d][ss][rn]
                File input_db = new File(ns.get("input_db").toString());
                File output_path = new File(ns.get("output_path").toString());
                File data_file = new File(ns.get("input_data_file").toString());
                String subsetSize = ns.getString("subset_size");
                String repetitionNum = ns.getString("repetition_num");
                String depth = ns.getString("depth");
                Integer d = 0;
                try{
                    d = Integer.valueOf(depth);
                }catch (NumberFormatException e){
                    System.out.println("Depth must be an integer");
                }
                Integer ss = Integer.valueOf(subsetSize);
                Integer rn = Integer.valueOf(repetitionNum);
                MeasuredDatabase mdb = new MeasuredDatabase(input_db, output_path);
                mdb.mouseEmpiricalNullDistribution(data_file,d, ss, rn);
            }
            /////////////////////////////////////////////////////////////////////////////////
            else{
                System.out.println("Mode \""+ mode+"\" not recognized, Options are:\n" +
                        "                    \"CreateDB\", takes an OWL file [-iof], an output path [-op], an update boolean [-u] (can be T or F), and the species of graph you'd like to make [-s] (can be human (h) or mouse(m))\n" +
                        "                    \"PrintDatabase\", takes in input database [-idb]\n" +
                        "                    \"WriteAllUIDs\", takes an input database [-idb] and an output path [-op]\n" +
                        "                    \"WritePhos\", takes an input database [-idb] and an output path [-op]\n" +
                        "                    \"WriteDBtoSIF\", takes an input database [-idb] and an output path [-op]\n" +
                        "                    \"MapPeptides\", takes an input database [-idb] and an output path [-op], and a file to map onto the database [-idf]\n" +
                        "                    \"IntegratePSP\", takes in input Reactome database [-idb], takes the PSP database [-psp], and an output path [-op]\n" +
                        "                    \"BinomialNeighbourhoodAnalysis\", takes in a measured input database [-idb], an output path [-op], and the depth of the traversal [-d]\n"+
                        "                    \"ShortestPath\", takes in a measured input database [-idb], an output path [-op], a starting node id [-sid], a ending node id [-eid], and the experiment name of interest [-en]\n"+
                        "                    \"MinimalConnectionNetwork\", takes in a measured input database [-idb], an output path [-op]\n"+
                        "                    \"TraversalAnalysis\", takes in a measured input database [-idb], an output path [-op], a UniProt ID or database ID to look downstream of [-p], the direction of the traversal [-dir], and the experiment name of interest [-en]\n"

                );

            }
        }catch (ArgumentParserException e){
            parser.handleError(e);
            System.exit(1);
        }





    }

}
