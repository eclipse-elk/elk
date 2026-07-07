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

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Queue;

import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * Main class of the comment attachment framework. Comment attachment infers the relation between comments and the
 * objects they refer to.
 * 
 * <p>
 * The comment attachment framework mainly consists of this class and the following interfaces to make
 * the framework customizable to each use case:
 * </p>
 * <ul>
 *   <li>
 *     {@link IDataProvider}<br/>
 *     This is what gives the framework access to comments and the possible attachment targets they can be attached to.
 *     The interface supports hierarchical data structures (such as ELK graphs) and also provides the means to take
 *     appropriate actions upon finding that a comment and an attachment target are in fact to be attached.
 *   </li>
 *   <li>
 *     {@link IExplicitAttachmentProvider}<br/>
 *     Some visual languages allow developers to establish explicit attachments between comments and targets. An
 *     {@link IExplicitAttachmentProvider} knows how to retrieve such explicit attachments. Clients can configure
 *     whether the presence of explicit attachments in a data set should disable heuristic attachments. After all, if
 *     the developer already defined some explicit attachments, chances are they wanted the other comments to be
 *     unattached.
 *   </li>
 *   <li>
 *     {@link IBoundsProvider}<br/>
 *     Comment attachment often includes looking at the size and position of comments and attachment targets. A bounds
 *     provider knows how to provide these information.
 *   </li>
 *   <li>
 *     {@link IFilter}<br/>
 *     Some comments clearly refer to a diagram as a whole. Filters filter these out to keep them from being attached
 *     to anything. Examples are comments that identify the authors of a diagram.
 *   </li>
 *   <li>
 *     {@link IMatcher}<br/>
 *     Matchers provide a heuristic assessment as to whether a given comment-target pair is likely to be attached or
 *     not. Matchers can be based on a lot of different metrics, but have to provide a way of normalizing their results
 *     to {@code [0, 1]}. Attachment decisions are based on the values computed by matchers.
 *   </li>
 *   <li>
 *     {@link IDecider}<br/>
 *     Attachment deciders have the final say on which attachment target a comment will be attached to based on the
 *     information computed by matchers. Filtered comments will never reach a decider.
 *   </li>
 * </ul>
 * 
 * <h3>Configuration</h3>
 * <p>
 * All of the different aspects of comment attachment can be configured through the configuration methods provided by
 * this class. Each configuration method documents the default used if it is not called. The following configuration
 * methods will probably be of particular interest:
 * </p>
 * <ul>
 *   <li>
 *     {@link #withExplicitAttachmentProvider(IExplicitAttachmentProvider)}<br/>
 *     Even if no heuristics are configured, explicit attachments may be of interest.
 *   </li>
 *   <li>
 *     {@link #addFilter(IFilter)}<br/>
 *     Filtering out standalone comments will probably be of interest.
 *   </li>
 *   <li>
 *     {@link #addHeuristic(IMatcher)}<br/>
 *     Heuristics are the backbone of the comment attachment framework. Unless you only want to add explicit
 *     attachments, you will want to add heuristics to base attachment decisions on.
 *   </li>
 *   <li>
 *     {@link #withBoundsProvider(IBoundsProvider)}<br/>
 *     If a single bounds provider instance is used in the whole system and that needs to preprocess things, register
 *     it using this method to have the comment attacher call its preprocessing method at the appropriate time.
 *   </li>
 * </ul>
 * 
 * @param <C>
 *            type of comments.
 * @param <T>
 *            type of attachment targets.
 */
public final class CommentAttacher<C, T> {
    
    /** Whether to go down on the graph or not. */
    private boolean includeHierarchy = true;
    /** Whether the presence of explicit attachments disables the heuristic attachments. */
    private boolean explicitAttachmentsDisableHeuristics = true;
    /** Retrieves graph elements explicitly attached to a comment by the user. */
    private IExplicitAttachmentProvider<C, T> explicitAttachmentProvider = (a) -> null;
    /** The bounds provider to be used. */
    private IBoundsProvider<C, T> boundsProvider = null;
    /** List of filters. */
    private List<IFilter<C>> filters = Lists.newArrayList();
    /** List of matchers. */
    private List<IMatcher<C, T>> matchers = Lists.newArrayList();
    /** The attachment decider. */
    private IDecider<T> decider = new AggregatedMatchDecider<>();
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    /**
     * Configures whether comment attachment should only attach comments on the current hierarchy level.
     * 
     * <p>
     * If this method is not called, comments on all hierarchy levels below the one passed to the attacher are attached.
     * </p>
     * 
     * @param limit
     *            {@code true} if comment attachment should be limited to the current hierarchy level.
     * @return this comment attacher for method chaining.
     */
    public CommentAttacher<C, T> limitToCurrentHierarchyLevel(final boolean limit) {
        includeHierarchy = !limit;
        return this;
    }
    
    /**
     * Configures whether comment attachment should keep using matchers to attach comments even if explicit attachments
     * are found.
     * 
     * <p>
     * If this method is not called, heuristics are disabled if explicit attachments are found.
     * </p>
     * 
     * @param keepEnabled
     *            {@code true} if heuristically found attachments should be applied even when explicit attachments are
     *            found.
     * @return this comment attacher for method chaining.
     */
    public CommentAttacher<C, T> keepHeuristicsEnabledWithExplicitAttachments(final boolean keepEnabled) {
        explicitAttachmentsDisableHeuristics = !keepEnabled;
        return this;
    }
    
    /**
     * Configures comment attachment to use the given explicit attachment provider. See the class documentation for more
     * information on how explicit attachments are handled.
     * 
     * <p>
     * If this method is not called, no explicit attachments are recognized by default.
     * </p>
     * 
     * @param provider
     *            the explicit attachment provider. Can be {@code null} to disable explicit attachments.
     * @return this comment attacher for method chaining.
     */
    public CommentAttacher<C, T> withExplicitAttachmentProvider(final IExplicitAttachmentProvider<C, T> provider) {
        if (provider == null) {
            explicitAttachmentProvider = (a) -> null;
        } else {
            explicitAttachmentProvider = provider;
        }
        
        return this;
    }
    
    /**
     * Configures comment attachment to call {@link IBoundsProvider#preprocess(IDataProvider, boolean)) and
     * {@link IBoundsProvider#cleanup()} on the given bounds provider at the appropriate times.
     * 
     * <p>
     * If this method is not called, no bounds provider will be notified of preprocessing and cleanup.
     * </p>
     * 
     * @param provider
     *            the non-{@code null} bounds provider.
     * @return this comment attacher for method chaining.
     * @throws IllegalArgumentException
     *             if {@code provider == null}.
     */
    public CommentAttacher<C, T> withBoundsProvider(final IBoundsProvider<C, T> provider) {
        Objects.requireNonNull(provider, "The bounds provider must not be null.");

        boundsProvider = provider;
        return this;
    }
    
    /**
     * Adds the given filter to the list of filters to be used. See the class documentation for more information on
     * filters.
     * 
     * <p>
     * If this method is not called, all comments are considered eligible for attachment.
     * </p>
     * 
     * @param filter
     *            the non-{@code null} filter.
     * @return this comment attacher for method chaining.
     * @throws IllegalArgumentException
     *             if {@code filter == null}.
     */
    public CommentAttacher<C, T> addFilter(final IFilter<C> filter) {
        Objects.requireNonNull(filter, "The filter must not be null.");

        filters.add(filter);
        return this;
    }
    
    /**
     * Adds the given matcher to the list of matchers to be used. See the class documentation for more information on
     * matchers.
     * 
     * <p>
     * If this method is not called, no comment will be heuristically attached to anything.
     * </p>
     * 
     * @param matcher
     *            the non-{@code null} matcher.
     * @return this comment attacher for method chaining.
     * @throws IllegalArgumentException
     *             if {@code heuristic == null}.
     */
    public CommentAttacher<C, T> addMatcher(final IMatcher<C, T> matcher) {
        Objects.requireNonNull(matcher, "The matcher must not be null.");

        matchers.add(matcher);
        return this;
    }
    
    /**
     * Configures comment attachment to use the given decider. See the class documentation for more information on
     * deciders.
     * 
     * <p>
     * If this method is not called, no comments are ever heuristically attached to anything. This is probably not what
     * you want.
     * </p>
     * 
     * @param attachmentDecider
     *            the non-{@code null} decider.
     * @return this comment attacher for method chaining.
     * @throws IllegalArgumentException
     *             if {@code decider == null}.
     */
    public CommentAttacher<C, T> withAttachmentDecider(final IDecider<T> attachmentDecider) {
        Objects.requireNonNull(attachmentDecider, "The attachment decider must not be null.");

        this.decider = attachmentDecider;
        return this;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Actual Logic
    
    /**
     * Runs the comment attachment heuristic with its current configuration on the given data provider. Attaching
     * comments and targets is done by calling {@link IDataProvider#attach(Object, Object)}.
     * 
     * @param dataProvider
     *            provider for the data set to run comment attachment on.
     */
    public void attachComments(final IDataProvider<C, T> dataProvider) {
        // Do necessary preprocessing
        preprocess(dataProvider);
        
        // We keep track of all comment->target pairs we find (we may want to store the data provider that was used
        // to obtain the comment and the attachment target at some point)
        Collection<Pair<C, T>> explicitAttachments = Lists.newArrayList();
        Collection<Pair<C, T>> heuristicAttachments = Lists.newArrayList();
        
        // Iterate over all comments we should iterate over
        Queue<IDataProvider<C, T>> processingQueue = Lists.newLinkedList();
        processingQueue.add(dataProvider);
        
        while (!processingQueue.isEmpty()) {
            IDataProvider<C, T> currentHierarchy = processingQueue.poll();
            
            for (C comment : currentHierarchy.provideComments()) {
                // Find explicit attachment
                T explicitAttachment = explicitAttachmentProvider.findExplicitAttachment(comment);
                
                if (explicitAttachment != null) {
                    // We'll use the explicit attachment
                    explicitAttachments.add(Pair.of(comment, explicitAttachment));
                    
                } else if ((explicitAttachments.isEmpty() || !explicitAttachmentsDisableHeuristics)
                        && isEligibleForHeuristicAttachment(comment)) {
                    
                    // We don't even bother calculating the heuristics if they are disabled by explicit
                    // attachments anyway. If that's not the case, let's roll!
                    T heuristicAttachment = findMatch(currentHierarchy, comment);
                    
                    if (heuristicAttachment != null) {
                        heuristicAttachments.add(Pair.of(comment, heuristicAttachment));
                    }
                }
            }
            
            // If we include hierarchy, check if the current node has children
            if (includeHierarchy) {
                processingQueue.addAll(currentHierarchy.provideSubHierarchies());
            }
        }
        
        // Tell everyone to clean up after themselves
        cleanup();
        
        // Attach what we've found
        edgeifyFoundAttachments(dataProvider, explicitAttachments, heuristicAttachments);
    }

    /**
     * Gives everything a chance to preprocess stuff.
     * 
     * @param dataProvider
     *            provider for the data set to run comment attachment on.
     */
    private void preprocess(final IDataProvider<C, T> dataProvider) {
        explicitAttachmentProvider.preprocess(dataProvider, includeHierarchy);
        filters.stream().forEach((f) -> f.preprocess(dataProvider, includeHierarchy));
        matchers.stream().forEach((h) -> h.preprocess(dataProvider, includeHierarchy));
        
        if (boundsProvider != null) {
            boundsProvider.preprocess(dataProvider, includeHierarchy);
        }
    }

    /**
     * Checks whether the given comment is eligible for heuristically finding attachment targets.
     * 
     * @param comment
     *            the comment to check.
     * @return {@code true} if the comment may be attached to things.
     */
    private boolean isEligibleForHeuristicAttachment(final C comment) {
        return filters.stream().allMatch((f) -> f.eligibleForAttachment(comment));
    }
    
    /**
     * Runs the matxhers on the given comment and returns the target it is to be attached to, if any.
     * 
     * @param dataProvider
     *            provider to ask for possible attachment targets.
     * @param comment
     *            the comment to run the matchers on.
     * @return the taget the comment is to be attached to, or {@code null} if there isn't any.
     */
    @SuppressWarnings("unchecked")
    private T findMatch(final IDataProvider<C, T> dataProvider, final C comment) {
        // If there are no heuristics, return nothing...
        if (matchers.isEmpty()) {
            return null;
        }
        
        // Collect attachment target candidates
        Collection<T> candidates = dataProvider.provideTargetsFor(comment);
        if (candidates == null || candidates.isEmpty()) {
            return null;
        }
        
        // Collect the matcher results in this map, indexed by attachment target, then indexed by the matcher
        Map<T, Map<Class<? extends IMatcher<?, T>>, Double>> results = Maps.newHashMap();
        
        for (T candidate : candidates) {
            Map<Class<? extends IMatcher<?, T>>, Double> candidateResults = Maps.newHashMap();
            results.put(candidate, candidateResults);
            
            // Run the normalized heuristics and collect their results in an array
            for (IMatcher<C, T> heuristic : matchers) {
                candidateResults.put(
                        (Class<? extends IMatcher<?, T>>) heuristic.getClass(),
                        heuristic.normalized(comment, candidate));
            }
        }
        
        // Decide which attachment target to attach the comment to
        return decider.makeAttachmentDecision(results);
    }
    
    /**
     * Tells the data provider to attach comments to their attachment targets. If there are explicit attachments, the
     * heuristic attachments are only applied if the presence of explicit attachments doesn't disable the heuristics.
     * 
     * @param dataProvider
     *            provider to delegate the actual attachment process to.
     * @param explicitAttachments
     *            collection of explicit attachments from comment nodes to graph elements.
     * @param heuristicAttachments
     *            collection of heuristic attachments from comment nodes to graph elements.
     */
    private void edgeifyFoundAttachments(final IDataProvider<C, T> dataProvider,
            final Collection<Pair<C, T>> explicitAttachments, final Collection<Pair<C, T>> heuristicAttachments) {
        
        // Explicit attachments
        for (Pair<C, T> explicitAttachment : explicitAttachments) {
            dataProvider.attach(explicitAttachment.getFirst(), explicitAttachment.getSecond());
        }
        
        // Only apply heuristic attachments if they are still enabled
        if (explicitAttachments.isEmpty() || !explicitAttachmentsDisableHeuristics) {
            for (Pair<C, T> heuristicAttachment : heuristicAttachments) {
                dataProvider.attach(heuristicAttachment.getFirst(), heuristicAttachment.getSecond());
            }
        }
    }
    
    /**
     * Gives everything a chance to cleanup.
     */
    private void cleanup() {
        explicitAttachmentProvider.cleanup();
        filters.stream().forEach((f) -> f.cleanup());
        matchers.stream().forEach((h) -> h.cleanup());
        
        if (boundsProvider != null) {
            boundsProvider.cleanup();
        }
    }
    
}
