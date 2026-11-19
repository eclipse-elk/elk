/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.transform;

/**
 * Definition of Neato distance computation models.
 *
 * @author msp
 */
public enum NeatoModel {
    
    /** Use length of the shortest path. */
    SHORTPATH,
    /** Use the circuit resistance model. */
    CIRCUIT,
    /** Use the subset model. */
    SUBSET;
    
    /**
     * Parse the given string into a Neato model, ignoring case.
     * 
     * @param string a string
     * @return the corresponding model
     */
    public static NeatoModel parse(final String string) {
        return valueOf(string.toUpperCase());
    }
    
    /**
     * Return the literal value as understood by Graphviz.
     * 
     * @return the literal value
     */
    public String literal() {
        return super.toString().toLowerCase();
    }
    
}
