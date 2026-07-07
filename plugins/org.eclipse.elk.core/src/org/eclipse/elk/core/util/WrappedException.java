/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

/**
 * A runtime exception that can be used to wrap checked exceptions. Use this where it is
 * appropriate to forward an error to the next point where it can be handled (i.e. displayed
 * to the user) without the need to explicitly declare the error in every method signature.
 *
 * @author msp
 */
public class WrappedException extends RuntimeException {

    /** the serial version UID. */
    private static final long serialVersionUID = -1630132187697677735L;
    
    /**
     * Create a wrapped exception.
     * 
     * @param cause the error that caused this exception
     */
    public WrappedException(final Throwable cause) {
        super(cause);
    }
    
    /**
     * Create a wrapped exception with additional message.
     * 
     * @param message an additional message for information
     * @param cause the error that caused this exception
     */
    public WrappedException(final String message, final Throwable cause) {
        super(message, cause);
    }

}
