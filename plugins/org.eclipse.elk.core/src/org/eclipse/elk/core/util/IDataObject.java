/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.io.Serializable;

/**
 * An interface for data types, which should be serializable using {@link #toString()} and
 * parsable using {@link #parse(String)}. The default constructor must always be
 * accessible and create an instance with default content.
 *
 * @author msp
 */
public interface IDataObject extends Serializable {
    
    /**
     * Parse the given string and set the content of this data object.
     * 
     * @param string a string
     * @throws IllegalArgumentException if the string does not have the expected format
     */
    void parse(String string);

}
