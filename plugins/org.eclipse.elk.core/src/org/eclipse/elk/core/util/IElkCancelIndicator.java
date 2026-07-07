/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

/**
 * Interface for indicators of task cancelation. Possibly long-running tasks such as layout algorithms
 * should query the indicator from time to time in order to allow users to abort the process.
 */
public interface IElkCancelIndicator {

    /**
     * Returns whether cancellation of the task has been requested.
     * 
     * @return true if cancellation has been requested
     */
    boolean isCanceled();

}
