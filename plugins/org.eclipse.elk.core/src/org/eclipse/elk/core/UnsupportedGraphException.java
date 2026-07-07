/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core;

/**
 * Thrown when a layout algorithm is executed on a graph that is not supported.
 *
 * @author msp
 */
public class UnsupportedGraphException extends RuntimeException {

    /** the serial version UID. */
    private static final long serialVersionUID = 669762537737088914L;
    
    /**
     * Create an unsupported graph exception with no parameters.
     */
    public UnsupportedGraphException() {
        super();
    }
    
    /**
     * Create an unsupported graph exception with a message.
     * 
     * @param message a message
     */
    public UnsupportedGraphException(final String message) {
        super(message);
    }
    
    /**
     * Create an unsupported graph exception with a message and a cause.
     * 
     * @param message a message
     * @param cause a cause
     */
    public UnsupportedGraphException(final String message, final Throwable cause) {
        super(message, cause);
    }

}
