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
 * A matcher determines how likely it is that a comment is related to an attachment target. Matchers provide a raw
 * heuristic result and a normalized result. The latter is constrained to be in {@code [0, 1]} (including the
 * boundaries), while the former is unconstrained, but may choose to follow the same constraints.
 * 
 * <p>
 * Note that if an implementation holds state, all resources should be released once {@link #cleanup()} is called.
 * </p>
 * 
 * @param <C>
 *            type of comments.
 * @param <T>
 *            type of attachment targets.
 */
public interface IMatcher<C, T> {
    
    /**
     * Computes the raw heuristic value for the given comment and the given attachment target. The heuristic value is
     * expressed as an arbitrary number.
     * 
     * <p>
     * Usually, the overall comment attachment will be based on the normalized value given by
     * {@link #normalized(Object, Object)} instead of on the raw value, but the raw value may be interesting for
     * evaluation purposes. However, the raw value may well be equal to the normalized value. An example are heuristics
     * that conceptually return boolean values: a pair of attachment targets either is considered to be related
     * ({@code 1}) or not ({@code 0}.
     * </p>
     * 
     * @param comment
     *            the comment node.
     * @param target
     *            the attachment target.
     * @return the raw heuristic result for the given comment and graph element.
     */
    double raw(C comment, T target);
    
    /**
     * Computes the normalized heuristic value in {@code [0, 1]} (including the boundaries) for the given comment and
     * the given attachment target.
     * 
     * @param comment
     *            the comment node.
     * @param target
     *            the attachment target.
     * @return the normalized heuristic result for the given comment and attachment target.
     */
    double normalized(C comment, T target);
    
    /**
     * Does any preprocessing necessary. This method is called before the first invocation of the non-default interface
     * methods for a given graph.
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
     * called after the last invocation of the non-default interface methods for a given ata set.
     * 
     * @implSpec The default implementation does nothing.
     */
    default void cleanup() {
    }
    
}
