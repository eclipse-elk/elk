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

/**
 * Knows which attachment targets comments were explicitly attached to by the user. This of course only makes sense for
 * comment attachment in languages whose editors allow users to explicitly attach comments to things.
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
 * @param <T>
 *            type of attachment targets.
 */
@FunctionalInterface
public interface IExplicitAttachmentProvider<C, T> {

    /**
     * Finds and returns the attachment target the given comment was explicitly attached to by the user, if any.
     * 
     * @param comment
     *            the comment whose explicit attachment to find.
     * @return the target the comment is explicitly attached to, or {@code null} if there is none. If the target is not
     *         {@code null}, it must be in the same hierarchy level as the comment node.
     */
    T findExplicitAttachment(C comment);

    /**
     * Does any preprocessing necessary. This method is called before the first invocation of
     * {@link #findExplicitAttachment(Object)} for a given graph.
     * 
     * @implSpec
     * The default implementation does nothing.
     * 
     * @param dataProvider
     *            the data set provider.
     * @param includeHierarchy
     *            if {@code true}, the preprocessing is done not only on the current hierarchy level, but also on all
     *            sub levels. Implementations may choose to behave differently depending on this value.
     */
    default void preprocess(IDataProvider<C, T> dataProvider, boolean includeHierarchy) {
    }

    /**
     * Does any cleaning necessary to get the implementation ready for the next comment attachment run. This method is
     * called after the last invocation of {@link #findExplicitAttachment(Object)} for a given data set.
     * 
     * @implSpec The default implementation does nothing.
     */
    default void cleanup() {
    }

}
