/*******************************************************************************
 * Copyright (c) 2010, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.graph;


/**
 * A label in the layered graph structure.
 */
public final class LLabel extends LShape {
    
    /** the serial version UID. */
    private static final long serialVersionUID = -264988654527750053L;
    
    /** text of the label. */
    private String text;
    
    /**
     * Creates a label with empty text.
     */
    public LLabel() {
        this("");
    }
    
    /**
     * Creates a label.
     * 
     * @param thetext text of the label
     */
    public LLabel(final String thetext) {
        this.text = thetext;
    }
    
    /**
     * Returns the text of the label.
     * 
     * @return the text
     */
    public String getText() {
        return text;
    }

    @Override
    public String toString() {
        String designation = getDesignation();
        if (designation == null) {
            return "label";
        } else {
            return "l_" + designation;
        }
    }

    @Override
    public String getDesignation() {
        if (!(text == null || text.isEmpty())) {
            return text;
        }
        return super.getDesignation();
    }
    
}
