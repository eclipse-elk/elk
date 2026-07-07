/*******************************************************************************
 * Copyright (c) 2012, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.properties;

/**
 * Graph features used for automatic recognition of the suitability of layout algorithms.
 * 
 * <p><em>Note:</em>
 * Originally, this enumeration was to be found in the {@code properties} package of the core plug-in. However,
 * the meta data compiler needs access to this enumeration. Making it depend on the core plug-in introduces an
 * unfortunate cyclic dependency: to compile the meta-data compiler, we need to compile the core plug-in. To
 * compile the core plug-in, we need the meta-data compiler.</p>
 *
 * @author msp
 * @author cds
 */
public enum GraphFeature {

    /** Edges connecting a node with itself. */
    SELF_LOOPS("Edges connecting a node with itself."),
    /** Self-loops routed through a node instead of around it. */
    INSIDE_SELF_LOOPS("Self-loops routed through a node instead of around it."),
    /** Multiple edges with the same source and target node. */
    MULTI_EDGES("Multiple edges with the same source and target node."),
    /** Labels that are associated with edges. */
    EDGE_LABELS("Labels that are associated with edges."),
    /** Edges are connected to nodes over ports. */
    PORTS("Edges are connected to nodes over ports."),
    /**
     * Edges that connect nodes from different hierarchy levels and are incident to compound nodes.
     * @see LayoutOptions#LAYOUT_HIERARCHY
     */
    COMPOUND("Edges that connect nodes from different hierarchy levels"
            + "and are incident to compound nodes."),
    /** Edges that connect nodes from different clusters, but not the cluster parent nodes. */
    CLUSTERS("Edges that connect nodes from different clusters, but not the cluster parent nodes."),
    /** Multiple connected components. */
    DISCONNECTED("Multiple connected components.");
    
    /** The description of the graph feature. */
    private String description;
    
    /**
     * Create a new graph feature with no description.
     */
    GraphFeature() {
        description = "";
    }
    
    /**
     * Creates a new graph feature with the given description.
     * 
     * @param description
     *            the description of the graph feature.
     */
    GraphFeature(final String description) {
        this.description = description;
    }

    /**
     * Returns the description of the graph feature.
     * 
     * @return the description
     */
    public String getDescription() {
        return description;
    }
    

}
