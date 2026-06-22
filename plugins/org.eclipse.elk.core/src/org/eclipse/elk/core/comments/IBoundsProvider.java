/*******************************************************************************
 * Copyright (c) 2016, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.comments;

import java.awt.geom.Rectangle2D;
import java.awt.geom.Rectangle2D.Double;
import java.util.HashMap;
import java.util.Map;

/**
 * Knows how to retrieve the position and size of comments and attachment targets.
 * 
 * <p>
 * Note that if an implementation holds state, all resources should be released once {@link #cleanup()} is called.
 * </p>
 * 
 * <p>
 * Bounds providers can be turned into caching bounds providers by calling {@link #cached()} on a bounds provider.
 * Cached bounds providers only calculate the bounds of nodes the first time, and use cached bounds afterwards. Use
 * cached bounds providers as proxies for computationally expensive bounds providers. It is usually a good idea to use
 * a single cached bounds provider for the whole comment attachment system.
 * </p>
 * 
 * @param <C>
 *            type of comments.
 * @param <T>
 *            type of attachment targets.
 * @see ElkGraphBoundsProvider
 */
public interface IBoundsProvider<C, T> {

    /**
     * Returns the position and size to be used for the given attachment target during comment attachment.
     * 
     * @param comment
     *            the comment whose bounds to retrieve.
     * @return the comment's bounds, or {@code null} if the provider was unable to retrieve them.
     */
    Rectangle2D.Double boundsForComment(C comment);
    
    /**
     * Returns the position and size to be used for the given attachment target during comment attachment.
     * 
     * @param target
     *            the target whose bounds to retrieve.
     * @return the target's bounds, or {@code null} if the provider was unable to retrieve them.
     */
    Rectangle2D.Double boundsForTarget(T target);
    
    /**
     * Does any preprocessing necessary. This method is called before the first invocation of the bounds methods for a
     * given data set.
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
     * called when no more invocations of {@link #boundsForComment(Object)} and {@link #boundsForTarget(Object)} will
     * follow for a given data set.
     * 
     * @implSpec The default implementation does nothing.
     */
    default void cleanup() {
    }
    
    /**
     * Returns a bounds provider that acts as a cache around this bounds provider.
     * 
     * @return a caching bounds provider around this bounds provider.
     */
    default IBoundsProvider<C, T> cached() {
        return new IBoundsProvider<C, T>() {

            /** Cache for comment bounds. */
            private final Map<C, Rectangle2D.Double> commentBoundsCache = new HashMap<>();
            /** Cache for target bounds. */
            private final Map<T, Rectangle2D.Double> targetBoundsCache = new HashMap<>();
            
            @Override
            public Double boundsForComment(final C comment) {
                // The cache may actually map a comment to null if the bounds provider was unable to
                // retrieve bounds for the comment
                if (commentBoundsCache.containsKey(comment)) {
                    return commentBoundsCache.get(comment);
                } else {
                    Rectangle2D.Double bounds = IBoundsProvider.this.boundsForComment(comment);
                    commentBoundsCache.put(comment, bounds);
                    return bounds;
                }
            }
            
            @Override
            public Double boundsForTarget(final T target) {
                // The cache may actually map a target to null if the bounds provider was unable to
                // retrieve bounds for the target
                if (targetBoundsCache.containsKey(target)) {
                    return targetBoundsCache.get(target);
                } else {
                    Rectangle2D.Double bounds = IBoundsProvider.this.boundsForTarget(target);
                    targetBoundsCache.put(target, bounds);
                    return bounds;
                }
            }
            
            @Override
            public void preprocess(final IDataProvider<C, T> dataProvider, final boolean includeHierarchy) {
                // Forward to wrapped provider
                IBoundsProvider.this.preprocess(dataProvider, includeHierarchy);
            }

            @Override
            public void cleanup() {
                commentBoundsCache.clear();
                targetBoundsCache.clear();
                
                // Forward to wrapped provider
                IBoundsProvider.this.cleanup();
            }
            
        };
    }
    
}
