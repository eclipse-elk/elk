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
 * Definition of available Graphviz commands.
 *
 * @author msp
 */
public enum Command {
    
    /** invalid command. */
    INVALID,
    /** command for Dot layout. */
    DOT,
    /** command for Neato layout. */
    NEATO,
    /** command for Twopi layout. */
    TWOPI,
    /** command for Fdp layout. */
    FDP,
    /** command for Circo layout. */
    CIRCO;

    /**
     * Parse the given string into a command, ignoring case.
     * 
     * @param string a string
     * @return the corresponding command, or {@code INVALID} if the string
     *     does not hold a valid command
     */
    public static Command parse(final String string) {
        try {
            return valueOf(string.toUpperCase());
        } catch (IllegalArgumentException exception) {
            return INVALID;
        }
    }
    
    @Override
    public String toString() {
        return super.toString().toLowerCase();
    }
    
}
