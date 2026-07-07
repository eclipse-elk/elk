/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.data;

/**
 * Interface for data holder classes that describe the layout infrastructure.
 * Implementing classes are used to store meta information and are accessed using
 * {@link LayoutMetaDataService}.
 *
 * @author msp
 */
public interface ILayoutMetaData {
    
    /**
     * Returns the identifier, that is a unique string used to look up the layout data instance.
     * 
     *  @return the identifier
     */
    String getId();
    
    /**
     * Returns the name, that is a user-friendly title that can be displayed in user interfaces.
     * While the name should generally be short, longer text can be written in the description.
     * In contrast to the identifier, the name can be ambiguous.
     *
     * @return the name
     */
    String getName();
    
    /**
     * Returns the description for display in user interfaces.
     *
     * @return the description
     */
    String getDescription();
    
}
