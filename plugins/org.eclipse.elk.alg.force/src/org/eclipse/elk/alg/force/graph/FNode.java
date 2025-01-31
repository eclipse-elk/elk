/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.force.graph;

import java.util.ArrayList;
import java.util.List;

/**
 * A node in the force graph.
 * 
 * @author owo
 * @author msp
 */
public final class FNode extends FParticle {

    /** the serial version UID. */
    private static final long serialVersionUID = 8663670492984978893L;

    // CHECKSTYLEOFF VisibilityModifier
    /** the identifier number. */
    public int id;
    // CHECKSTYLEON VisibilityModifier
    /** the node label. */
    private String label;
    /** The parent node. */
    private FNode parent;
    /** List of child nodes. */
    private List<FNode> children;
    
    /**
     * Create a new node.
     */
    public FNode() {
    }
    
    /**
     * Create a new node with given label.
     * 
     * @param label the label text
     */
    public FNode(final String label) {
        this.label = label;
    }

    /**
     * Create a new node with given parent node.
     * 
     * @param label the label text
     * @param theParent the parent node
     */
    public FNode(final String label, final FNode theParent) {
        this.label = label;
        this.parent = theParent;
    }
    
    @Override
    public String toString() {
        if (label == null || label.length() == 0) {
            return "n_" + id;
        } else {
            return "n_" + label;
        }
    }

    /**
     * Returns whether this node is a compound node.
     * 
     * @return true if compound
     */
    public boolean isCompound() {
        return children != null && children.size() > 0;
    }
    
    /**
     * Returns the label text of this node.
     * 
     * @return the label
     */
    public String getLabel() {
        return label;
    }

    /**
     * Returns the parent node.
     * 
     * @return the parent
     */
    public FNode getParent() {
        return parent;
    }

    /**
     * Returns the list of children, creating it if necessary.
     * 
     * @return the children
     */
    public List<FNode> getChildren() {
        if (children == null) {
            children = new ArrayList<>();
        }
        return children;
    }

    /**
     * Returns the depth of this node in the compound hierarchy.
     * 
     * @return the depth
     */
    public int getDepth() {
        int depth = 0;
        FNode node = parent;
        while (node != null) {
            node = node.getParent();
            depth++;
        }
        return depth;
    }
    
}
