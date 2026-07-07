/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.options;

/**
 * Enumeration of node alignment options. To be accessed using {@link CoreOptions#ALIGNMENT}.
 *
 * @author msp
 */
public enum Alignment {
    
    /** automatic alignment (default). */
    AUTOMATIC,
    /** left alignment. */
    LEFT,
    /** right alignment. */
    RIGHT,
    /** top alignment. */
    TOP,
    /** bottom alignment. */
    BOTTOM,
    /** center alignment. */
    CENTER;

}
