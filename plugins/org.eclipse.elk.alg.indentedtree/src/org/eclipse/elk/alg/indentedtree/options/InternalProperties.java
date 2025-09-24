/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.indentedtree.options;

import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

/**
 * Internal Property definitions for the Indented Tree Layout Algorithm
 */
public final class InternalProperties {

    
    /**
     * the total indentation of a node, as calculated through a depth-first search
     */
    public static final IProperty<Double> DEPTH = new Property<Double>("depth", 0.0);
    
    /**
     * the indentation of the edge to the left of a node, inherited from its "parent"
     */
    public static final IProperty<Double> EDGE_INDENTATION = new Property<Double>("edge_indentation", 0.0);
}
