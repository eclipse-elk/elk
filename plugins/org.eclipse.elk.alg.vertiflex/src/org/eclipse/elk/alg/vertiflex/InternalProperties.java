/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex;

import org.eclipse.elk.alg.vertiflex.p2relative.OutlineNode;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

/**
 * The internal properties of the tree layouter. Properties that are used internally
 * but do not need to be public are stored here.
 *
 */
public final class InternalProperties {
    
    /**
     * Private constructor to prevent initialization.
     */
    private InternalProperties() { }

    /**
     * Defines the left outline of this subtree.
     */
    public static final IProperty<OutlineNode> LEFT_OUTLINE = new Property<OutlineNode>("LEFT_OUTLINE");
    
    /**
     * Defines the right outline of this subtree.
     */
    public static final IProperty<OutlineNode> RIGHT_OUTLINE = new Property<OutlineNode>("RIGHT_OUTLINE");
    
    /**
     * Defines the maximum depth of outlines. This is the lowest point of the outline in the tree layout.
     */
    public static final IProperty<Double> OUTLINE_MAX_DEPTH = new Property<Double>("OUTLINE_MAX_DEPTH");
    
    /** Defines the canvas of the tree. */
    public static final IProperty<Double> MIN_X = new Property<Double>("MIN_X");
    
    /** Defines the canvas of the tree. */
    public static final IProperty<Double> MAX_X = new Property<Double>("MAX_X");
    
    /** Defines the canvas of the tree. */
    public static final IProperty<Double> MIN_Y = new Property<Double>("MIN_Y");
    
    /** Defines the canvas of the tree. */
    public static final IProperty<Double> MAX_Y = new Property<Double>("MAX_Y");
    
    /** The root node of the graph. */
    public static final IProperty<ElkNode> ROOT_NODE = new Property<ElkNode>("root");
    
    /**
     * Defines the bendpoint of an edge.  
     */
    public static final IProperty<Double> EDGE_BEND_HEIGHT = new Property<Double>("EDGE_BEND_HEIGHT");
    
    /**
     * Stores the model order of nodes. Smaller values come before larger values.
     */
    public static final IProperty<Integer> NODE_MODEL_ORDER = new Property<Integer>("Node Model Order");
    
}
