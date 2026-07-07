/*******************************************************************************
 * Copyright (c) 2016, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.comments;

import java.util.Map;

/**
 * An attachment decider has the final say on which attachment target to attach a comment to, if any. The decision is
 * based on the normalized heuristic values produced by the matchers for each possible target.
 * 
 * @param <T>
 *            type of attachment targets.
 */
@FunctionalInterface
public interface IDecider<T> {
    
    /**
     * Decides which target to attach a comment to, if any.
     * 
     * @param normalizedHeuristics
     *            maps possible attachment targets to a map from matchers to the normalized heuristic values they
     *            produced.
     * @return the selected attachment target, or {@code null} if the comment should be left unattached.
     */
    T makeAttachmentDecision(Map<T, Map<Class<? extends IMatcher<?, T>>, Double>> normalizedHeuristics);
    
}
