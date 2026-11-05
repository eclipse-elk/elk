/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Provides access to the objects that represent comments and the objects that represent possible attachment targets.
 * The data provided by data providers is what the rest of the system is going to be working on to determine
 * attachments. If the data structure comments and targets rely in can have sub-hierarchies, access to those is granted
 * through providing further data providers.
 * 
 * <p>
 * By default, all targets are returned as possible attachment targets for all comments. However, clients may choose to
 * restrict the set of possible targets for any given comment by implementing {@link #provideTargetsFor(Object)}.
 * </p>
 * 
 * <p>
 * To obtain a cached version of any data provider, simply call {@link #cached()}.
 * </p>
 * 
 * @param <C>
 *            type of comments.
 * @param <T>
 *            type of attachment targets.
 */
public interface IDataProvider<C, T> {
    
    /**
     * Provides the set of comments that the comment attachment code should run on.
     * 
     * @return iterable over all comments.
     */
    Collection<C> provideComments();
    
    /**
     * Provides the set of all possible attachment targets that comments might be attached to.
     * 
     * @return iterable over all possible attachment targets.
     */
    Collection<T> provideTargets();
    
    /**
     * Provides all attachment targets the given comment can be attached to.
     * 
     * @implSpec
     * The default implementation simply delegates to {@link #provideTargets()}.
     * 
     * @param comment
     *            the comment to provide valid attachment targets for.
     * @return the set of valid attachment targets for the given comment.
     */
    default Collection<T> provideTargetsFor(final C comment) {
        return provideTargets();
    }
    
    /**
     * Provides access to sub-levels in a hierarchy. Each level is represented by a data provider.
     * 
     * @return set of data providers for sub-levels in the hierarchy.
     */
    Collection<IDataProvider<C, T>> provideSubHierarchies();
    
    /**
     * Tells the data provider to do whatever it feels necessary to attach the given comment to the given attachment
     * target.
     * 
     * @param comment the comment to attach to the target.
     * @param target the target to attach the comment to.
     */
    void attach(C comment, T target);
    
    /**
     * Returns a data provider that acts as a cache around this provider.
     * 
     * @return a cached data provider around this data provider.
     */
    default IDataProvider<C, T> cached() {
        return new IDataProvider<C, T>() {
            
            /** Cache for comments. */
            private List<C> commentsCache = null;
            /** Cache for targets. */
            private List<T> targetsCache = null;
            /** Cache for comment-specific targets. */
            private final Map<C, Set<T>> commentTargetsCache = new HashMap<>();
            /** Cache for sub-level data providers. */
            private List<IDataProvider<C, T>> subProviderCache = null;

            @Override
            public Collection<C> provideComments() {
                if (commentsCache == null) {
                    commentsCache = new ArrayList<>();
                    commentsCache.addAll(IDataProvider.this.provideComments());
                }
                
                return commentsCache;
            }

            @Override
            public Collection<T> provideTargets() {
                if (targetsCache == null) {
                    targetsCache = new ArrayList<>();
                    targetsCache.addAll(IDataProvider.this.provideTargets());
                }
                
                return targetsCache;
            }

            @Override
            public Collection<T> provideTargetsFor(final C comment) {
                if (commentTargetsCache.containsKey(comment)) {
                    return commentTargetsCache.get(comment);
                } else {
                    Collection<T> targets = IDataProvider.this.provideTargetsFor(comment);
                    if (!targets.isEmpty()) {
                        commentTargetsCache.put(comment, new HashSet<>(targets));
                    }
                    return targets;
                }
            }

            @Override
            public Collection<IDataProvider<C, T>> provideSubHierarchies() {
                if (subProviderCache == null) {
                    subProviderCache = new ArrayList<>();
                    
                    // Add cached versions of the providers
                    IDataProvider.this.provideSubHierarchies().stream()
                        .map(provider -> provider.cached())
                        .forEach(provider -> subProviderCache.add(provider));
                }
                
                return subProviderCache;
            }
            
            @Override
            public void attach(final C comment, final T target) {
                IDataProvider.this.attach(comment, target);
            }
            
        };
    }
    
}
