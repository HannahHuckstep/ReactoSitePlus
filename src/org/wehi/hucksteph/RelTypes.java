package org.wehi.hucksteph;

import org.neo4j.graphdb.RelationshipType;

public enum RelTypes implements RelationshipType {
    INPUT,
    OUTPUT,
    CONTROLS,
    CATALYSIS,
    ID_BELONGS_TO,
    PHOSPHORYLATION,
    MODIFICATION,
    COMPONENT,
    PATHWAY_COMPONENT,
    SUB_PATHWAY,
    SMALL_MOL_EDGE
}
