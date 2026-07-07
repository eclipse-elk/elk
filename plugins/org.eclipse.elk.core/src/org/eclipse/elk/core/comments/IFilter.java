/*******************************************************************************
 * Copyright (c) 2016, 2017 Kiel University and others.
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
package org.eclipse.elk.core.comments;

/**
 * Determines if a given comment is eligible for attachment to an attachment target. There are types of comments that
 * can very well stand on their own, without any need for attachment. A filter identifies such comments.
 * 
 * <p>
 * Note that if an implementation holds state, all resources should be released once {@link #cleanup()} is called.
 * </p>
 * 
 * <p>
 * If clients implementing this interface don't need any preprocessing or cleanup, this interface can be used as a
 * functional interface.
 * </p>
 * 
 * @param <C>
 *            type of comments.
 */
@FunctionalInterface
public interface IFilter<C> {
    
    /**
     * Checks if the given comment can be attached to an attachment target.
     * 
     * @param comment
     *            the comment.
     * @return {@code true} if the comment can be attached to graph elements, {@code false} if the comment is meant to
     *         be standalone.
     */
    boolean eligibleForAttachment(C comment);
    
    /**
     * Does any preprocessing necessary. This method is called before the first invocation of
     * {@link #eligibleForAttachment(Object)} for a given graph.
     * 
     * @implSpec The default implementation does nothing.
     * 
     * @param dataProvider
     *            the data set provider.
     * @param includeHierarchy
     *            if {@code true}, the preprocessing is done not only on the current hierarchy level, but also on all
     *            sub levels. Implementations may choose to behave differently depending on this value.
     */
    default void preprocess(IDataProvider<C, ?> dataProvider, boolean includeHierarchy) {
    }
    
    /**
     * Does any cleaning necessary to get the implementation ready for the next comment attachment run. This method is
     * called after the last invocation of {@link #eligibleForAttachment(Object)} for a given data set.
     * 
     * @implSpec The default implementation does nothing.
     */
    default void cleanup() {
    }
    
}
