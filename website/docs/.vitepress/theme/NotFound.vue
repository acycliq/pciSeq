<script setup>
import { withBase } from 'vitepress'
</script>

<template>
  <div class="nf">
    <p class="nf-code">404</p>
    <h1 class="nf-title">Page not found</h1>
    <p class="nf-text">
      The page you are looking for does not exist, or it may have moved.
    </p>
    <a class="nf-home" :href="withBase('/')">Return home</a>
  </div>
</template>

<style scoped>
.nf {
  max-width: 520px;
  margin: 0 auto;
  padding: 12rem 1.5rem 18rem;
  text-align: center;
}

.nf-code {
  font-size: 0.8125rem;
  font-weight: 600;
  letter-spacing: 0.1em;
  color: var(--vp-c-brand-1);
  margin: 0;
}

.nf-title {
  font-size: 1.6rem;
  font-weight: 600;
  letter-spacing: -0.01em;
  margin: 0.5rem 0 0.75rem;
}

.nf-text {
  color: var(--vp-c-text-2);
  line-height: 1.6;
  margin: 0 0 2rem;
}

.nf-home {
  display: inline-block;
  border: 1px solid var(--vp-c-brand-1);
  border-radius: 8px;
  padding: 0.5rem 1.25rem;
  font-size: 0.875rem;
  font-weight: 500;
  color: var(--vp-c-brand-1);
  transition: background-color 0.2s, color 0.2s;
}

.nf-home:hover {
  background-color: var(--vp-c-brand-1);
  color: var(--vp-c-white);
}
</style>
