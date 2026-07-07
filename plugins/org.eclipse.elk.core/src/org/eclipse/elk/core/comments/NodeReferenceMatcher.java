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

import java.awt.geom.Rectangle2D;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.graph.ElkNode;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * A binary matcher that checks if target names are referenced in comment texts. This can either be the case or not,
 * hence the term <em>binary heuristic</em>. Indeed, the matcher will always output either 0 or 1. It will only output 1
 * if a target's name appears in a comment's text, and if no other target's name appears in it. The following
 * configuration methods have to be called before using this heuristic:
 * <ul>
 *   <li>{@link #withCommentTextProvider(Function)}</li>
 *   <li>{@link #withTargetNameProvider(Function)}</li>
 * </ul>
 * <p>
 * The heuristic can optionally be configured to only consider matches to be matches if the distance between a comment
 * and a matching target does not exceed a configurable threshold. In that case, the
 * {@link #withBoundsProvider(IBoundsProvider)} method has to be called as well.
 * </p>
 * 
 * <p>
 * The matcher can operate in two modes: a strict and a fuzzy mode.
 * </p>
 * 
 * 
 * <h3>Strict Mode</h3>
 * <p>
 * In strict mode, target names have to appear exactly as they are in the text of a comment, as a separate word.
 * </p>
 * 
 * 
 * <h3>Fuzzy Mode</h3>
 * <p>
 * In fuzzy mode, the matcher allows the comment text more freedom in mentioning target names. Apart from ignoring case,
 * as in strict mode, we here also allow arbitrary whitespace (including newlines) between two parts of the name. A part
 * is divided from another part either through spaces or through camel case.
 * </p>
 * 
 * @param <C>
 *            type of comments.
 * @param <T>
 *            type of attachment targets.
 */
public class NodeReferenceMatcher<C, T> implements IMatcher<C, T> {
    
    /** Function used to retrieve a comment's text. */
    private Function<C, String> commentTextFunction = null;
    /** Function used to retrieve a target's name. */
    private Function<T, String> targetNameFunction = null;
    /** The bounds provider to use. */
    private IBoundsProvider<C, T> boundsProvider = null;
    /** The maixmum distance attached comments may be away from each other. */
    private double maxDistance = -1;
    /** Whether to use fuzzy mode when looking for occurrences of a target's name in a comment's text. */
    private boolean fuzzy = false;
    /** The comment-target attachments we've found during preprocessing. */
    private Map<C, T> foundAttachments = Maps.newHashMap();
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    /**
     * Configures the matcher to use the given function to obtain the text for a comment.
     * 
     * <p>
     * If this method is not called, the matcher will throw an exception during preprocessing.
     * </p>
     * 
     * @param f
     *            the function to use.
     * @return this object for method chaining.
     */
    public NodeReferenceMatcher<C, T> withCommentTextProvider(final Function<C, String> f) {
        Objects.requireNonNull(f, "Comment text function cannot be null.");
        
        commentTextFunction = f;
        return this;
    }
    
    /**
     * Configures the matcher to use the given function to obtain the name of an attachment target.
     * 
     * <p>
     * If this method is not called, the matcher will throw an exception during preprocessing.
     * </p>
     * 
     * @param f
     *            the function to use.
     * @return this object for method chaining.
     */
    public NodeReferenceMatcher<C, T> withTargetNameProvider(final Function<T, String> f) {
        Objects.requireNonNull(f, "Node name function cannot be null.");
        
        targetNameFunction = f;
        return this;
    }
    
    /**
     * Configures the matcher to be fuzzy when looking for occurrences of target names in a comment's text. If fuzzy
     * mode is off, only exact occurrences will be found.
     * 
     * <p>
     * If this method is not called, the heuristic will default to strict mode.
     * </p>
     * 
     * @return this object for method chaining.
     */
    public NodeReferenceMatcher<C, T> withFuzzyMatching() {
        fuzzy = true;
        return this;
    }
    
    /**
     * Configures the matcher to use the given maximum attachment distance. A comment is only attached to a referenced
     * target if the distance between the two doesn't exceed this distance.
     * 
     * <p>
     * If this method is not called, the matcher does not impose a distance restriction.
     * </p>
     * 
     * @param distance
     *            the maximum possible distance. Negative values disable the distance restriction.
     * @return this object for method chaining.
     */
    public NodeReferenceMatcher<C, T> withMaximumAttachmentDistance(final double distance) {
        this.maxDistance = distance;
        return this;
    }
    
    /**
     * Configures the matcher to use the given bounds provider to determine the bounds of comments and targets.
     * 
     * <p>
     * If this method is not called, but a maximum attachment distance is configured, the filter will throw an
     * exception during preprocessing.
     * </p>
     * 
     * @param provider
     *            the bounds provider to use.
     * @return this object for method chaining.
     */
    public NodeReferenceMatcher<C, T> withBoundsProvider(final IBoundsProvider<C, T> provider) {
        Objects.requireNonNull(provider, "Bounds provider must not be null.");
        
        this.boundsProvider = provider;
        return this;
    }
    
    /**
     * Checks whether the current configuration is valid.
     * 
     * @throws IllegalStateException
     *             if the configuration is invalid.
     */
    private void checkConfiguration() {
        if (commentTextFunction == null) {
            throw new IllegalStateException("A comment text function is required.");
        }
        
        if (targetNameFunction == null) {
            throw new IllegalStateException("A node name function is required.");
        }
        
        if (maxDistance >= 0 && boundsProvider == null) {
            throw new IllegalStateException(
                    "A bounds provider must be installed if a maximum attachment distance is set.");
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // IMatcher

    @Override
    public void preprocess(final IDataProvider<C, T> dataProvider, final boolean includeHierarchy) {
        checkConfiguration();
        
        // Build up maps of comment texts
        List<Pair<C, String>> commentTexts = Lists.newArrayList();
        
        for (C comment : dataProvider.provideComments()) {
            String commentText = commentTextFunction.apply(comment);
            
            if (!Strings.isNullOrEmpty(commentText)) {
                commentTexts.add(Pair.of(comment, commentText));
            }
        }

        // Build up maps of target names
        List<Pair<T, String>> targetNames = Lists.newArrayList();
        for (T target : dataProvider.provideTargets()) {
            String targetName = targetNameFunction.apply(target);
            
            if (!Strings.isNullOrEmpty(targetName)) {
                targetNames.add(Pair.of(target, targetName));
            }
        }
        
        // Go find matches!
        goFindMatches(commentTexts, targetNames);
        
        // Recurse into sub levels
        if (includeHierarchy) {
            for (IDataProvider<C, T> subProvider : dataProvider.provideSubHierarchies()) {
                preprocess(subProvider, true);
            }
        }
    }

    @Override
    public void cleanup() {
        foundAttachments.clear();
    }

    @Override
    public double raw(final C comment, final T target) {
        // Simply lookup the comment-target pair in our map
        return foundAttachments.get(comment) == target ? 1 : 0;
    }

    @Override
    public double normalized(final C comment, final T target) {
        return raw(comment, target);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Matching
    
    /**
     * Runs off and finds matches between comments and attachment targets. Depending on whether fuzzy mode is on or off,
     * this method uses different kinds of regular expressions to find target names in comments. A target and a comment
     * match if the target's name is contained in the comment's text, and if no other target is mentioned in the
     * comment's text.
     * 
     * <p>
     * Matches are recorded in {@link #foundAttachments}.
     * </p>
     * 
     * @param commentTexts
     *            list of pairs of comments and their text. The text is expected to not be {@code null}.
     * @param targetNames
     *            list of pairs of targets and their names. The name is expected to not be {@code null}.
     */
    private void goFindMatches(final List<Pair<C, String>> commentTexts, final List<Pair<T, String>> targetNames) {
        // Produce regular expression patterns for all target names
        List<Pair<T, Pattern>> targetRegexps = Lists.newArrayListWithCapacity(targetNames.size());
        for (Pair<T, String> targetNamePair : targetNames) {
            Pattern regexp = fuzzy
                    ? fuzzyRegexpFor(targetNamePair.getSecond())
                    : strictRegexpFor(targetNamePair.getSecond());
            targetRegexps.add(Pair.of(targetNamePair.getFirst(), regexp));
        }
        
        // Check each pair of comment and target
        for (Pair<C, String> commentTextPair : commentTexts) {
            T foundTarget = null;
            
            for (Pair<T, Pattern> targetRegexpPair : targetRegexps) {
                Matcher matcher = targetRegexpPair.getSecond().matcher(commentTextPair.getSecond());
                
                if (matcher.find()) {
                    // We only want to establish associations if a target is the only one mentioned in a comment
                    if (foundTarget == null) {
                        foundTarget = targetRegexpPair.getFirst();
                    } else {
                        foundTarget = null;
                        break;
                    }
                }
            }
            
            // Record the match, if any and if the maximum attachment distance allows
            if (foundTarget != null) {
                if (maxDistance < 0) {
                    // There is no distance restriction
                    foundAttachments.put(commentTextPair.getFirst(), foundTarget);
                } else {
                    // Calculate the distance between the target and the comment
                    Rectangle2D.Double commentBounds = boundsProvider.boundsForComment(commentTextPair.getFirst());
                    Rectangle2D.Double nodeBounds = boundsProvider.boundsForTarget(foundTarget);
                    
                    if (DistanceMatcher.distance(commentBounds, nodeBounds) <= maxDistance) {
                        foundAttachments.put(commentTextPair.getFirst(), foundTarget);
                    }
                }
            }
        }
    }
    
    /**
     * Produces a fuzzy regular expression pattern for the given target name. The pattern matches the following
     * appearances of the target name:
     * <ul>
     *   <li>
     *     a space character in the target name can be represented by one or more whitespace and line break characters
     *     in the comment text.
     *   </li>
     *   <li>
     *     if the target name is camelCased, each upper-case character preceded by a lower-case character can be
     *     prefixed by one or more whitespace and line break characters in the comment text.
     *   </li>
     * </ul>
     * 
     * @param targetName
     *            the target name.
     * @return regular expression pattern for fuzzy containation.
     */
    private static Pattern fuzzyRegexpFor(final String targetName) {
        String trimmedTargetName = targetName.trim();
        StringBuffer regexp = new StringBuffer(targetName.length() * 2);
        StringBuffer currentSegment = new StringBuffer(targetName.length());
        
        for (int i = 0; i < trimmedTargetName.length(); i++) {
            char currC = trimmedTargetName.charAt(i);
            
            if (Character.isUpperCase(currC)) {
                // If the previous character was lower-case, end the current segment and insert whitespace placeholders
                if (i > 0 && Character.isLowerCase(trimmedTargetName.charAt(i - 1))) {
                    regexp.append(Pattern.quote(currentSegment.toString()));
                    currentSegment = new StringBuffer(targetName.length());
                    
                    regexp.append("[\\h\\v]*");
                }
                
                currentSegment.append(currC);
            } else if (Character.isWhitespace(currC)) {
                // The first of a series of whitespace characters must be replaced by whitespace placeholders in the
                // regular expression, and the current segment ends
                if (i > 0 && !Character.isWhitespace(trimmedTargetName.charAt(i - 1))) {
                    regexp.append(Pattern.quote(currentSegment.toString()));
                    currentSegment = new StringBuffer(targetName.length());
                    
                    regexp.append("[\\h\\v]*");
                }
            } else {
                // It's neither upper-case, nor whitespace, so just add it to the current segment
                currentSegment.append(currC);
            }
        }
        
        regexp.append(Pattern.quote(currentSegment.toString()));
        
        return Pattern.compile("\\b" + regexp.toString() + "\\b",
                Pattern.DOTALL | Pattern.CASE_INSENSITIVE);
    }
    
    /**
     * Produces a strict regular expression pattern for the given target name. The pattern matches the target name if
     * it's not part of a longer word.
     * 
     * @param targetName
     *            the target name.
     * @return regular expression pattern for strict containation.
     */
    private static Pattern strictRegexpFor(final String targetName) {
        return Pattern.compile("\\b" + Pattern.quote(targetName) + "\\b", Pattern.DOTALL);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors
    
    /**
     * Returns a map that maps comments to targets they were attached to by this heuristic. The returned map is
     * meaningful only if it is called between calls to {@link #preprocess(ElkNode, boolean)} and {@link #cleanup()}.
     * Comments that are not attached to anything don't appear in the map.
     * 
     * @return mapping of comments to targets.
     */
    public Map<C, T> getAttachments() {
        return foundAttachments;
    }

}
