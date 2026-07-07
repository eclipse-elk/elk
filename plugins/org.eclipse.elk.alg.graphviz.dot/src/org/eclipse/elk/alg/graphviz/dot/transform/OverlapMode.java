/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.transform;

/**
 * Definition of overlap removal modes.
 *
 * @author msp
 */
public enum OverlapMode {
    
    /** overlaps are retained ("true" in the Graphviz options). */
    NONE,
    /** overlaps are removed by uniformly scaling in x and y. */
    SCALE,
    /** x and y are separately scaled to remove overlaps. */
    SCALEXY,
    /** Prism, a proximity graph-based algorithm, is used to remove node overlaps. */
    PRISM,
    /** the layout will be scaled down as much as possible without introducing any overlaps. */
    COMPRESS;
    
    /**
     * Parse the given string into an overlap mode, ignoring case.
     * 
     * @param string a string
     * @return the corresponding mode
     */
    public static OverlapMode parse(final String string) {
        if (string.toLowerCase().equals("true")) {
            return NONE;
        }
        return valueOf(string.toUpperCase());
    }
    
    /**
     * Return the literal value as understood by Graphviz.
     * 
     * @return the literal value
     */
    public String literal() {
        if (this == NONE) {
            return "true";
        }
        return super.toString().toLowerCase();
    }

}
